#!/usr/bin/env -S octave --persist

pkg load openems;
pkg load csxcad;


if (!isempty(argv))
  Sim_Path = argv{1}
  analysis_only = 1;
else
  analysis_only = 0;
  Sim_Path = datestr(clock(), 'yyyy-mm-ss_HH:MM:SS')
endif




%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-6; % specify everything in um
MSL_length = 20e3;
MSL_width = 600;
substrate_thickness = 254;
substrate_epr = 3.66;


f_max = 2e9;


%% coil properties
separation = 10e3
## wire_radius = 100
num_loops = 10
loop_radius = 3e3
lead_height = 3 * loop_radius
coil_res = 100
coil_mesh_res = loop_radius / 10




%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD();
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PEC' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
resolution = c0/(f_max*sqrt(substrate_epr))/unit /50

increase_max = 1.3
mesh.x = SmoothMeshLines([-separation/2, separation/2], coil_mesh_res);
mesh.x = SmoothMeshLines( [-MSL_length, mesh.x, MSL_length], resolution);

mesh.y = SmoothMeshLines([-loop_radius, loop_radius], coil_mesh_res);
mesh.y = SmoothMeshLines( [-15*MSL_width, mesh.y, 15*MSL_width], resolution);
mesh.z = SmoothMeshLines([substrate_thickness + lead_height - loop_radius, substrate_thickness + lead_height + loop_radius], coil_mesh_res);
mesh.z = SmoothMeshLines( [linspace(0,substrate_thickness,5), mesh.z, 2*lead_height], resolution);
CSX = DefineRectGrid( CSX, unit, mesh );

%% substrate
CSX = AddMaterial( CSX, 'RO4350B' );
CSX = SetMaterialProperty( CSX, 'RO4350B', 'Epsilon', substrate_epr );
start = [mesh.x(1),   mesh.y(1),   0];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness];
CSX = AddBox( CSX, 'RO4350B', 0, start, stop );

%% MSL port
CSX = AddMetal( CSX, 'PEC' );
portstart = [ mesh.x(1), -MSL_width/2, substrate_thickness];
portstop  = [ -separation/2,  MSL_width/2, 0];
[CSX,port{1}] = AddMSLPort(CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'FeedShift', 10*resolution, 'MeasPlaneShift',  MSL_length/3);

portstart = [mesh.x(end), -MSL_width/2, substrate_thickness];
portstop  = [separation/2 ,  MSL_width/2, 0];
[CSX,port{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  MSL_length/3 );

%% wire coil
count = 0
p(1,1) = -separation/2;
p(2,1) = 0;
p(3,1) = substrate_thickness;

## p(1,2) = -separation/2;
## p(2,2) = 0;
## p(3,2) = substrate_thickness + lead_height;

dt = 1.0 / coil_res;

count = 1;
for n=0:num_loops-1
    for m=0:coil_res
      count = count + 1;
      p(1,count) = -separation/2 + n*(separation / num_loops) + m*(separation / num_loops / coil_res);
      p(2,count) = loop_radius * cos(2*pi*dt*m);
      p(3,count) = substrate_thickness + lead_height + loop_radius * sin(2*pi*dt*m);
    end
end
## count = count+1
## p(1,count) = separation/2;
## p(2,count) = 0;
## p(3,count) = substrate_thickness + lead_height;
count = count+1;
p(1,count) = separation/2;
p(2,count) = 0;
p(3,count) = substrate_thickness;
CSX = AddCurve(CSX, 'PEC', 0, p);

%% write/show/run the openEMS compatible xml-file
Sim_CSX = 'msl.xml';


if (analysis_only == 0)
  %% create empty simulation folder
  [status, message, messageid] = mkdir( Sim_Path );
  WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
  CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
  RunOpenEMS( Sim_Path, Sim_CSX );
endif
 
%% post-processing
close all
f = linspace( 1e6, f_max, 1601 );
port = calcPort( port, Sim_Path, f, 'RefImpedance', 50);

s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;

plot(f/1e9,20*log10(abs(s11)),'LineWidth',2);
hold on;
grid on;
plot(f/1e9,20*log10(abs(s21)),'LineWidth',2);
legend('S_{11}','S_{21}');
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (GHz) \rightarrow','FontSize',12);
ylim([-40 2]);
