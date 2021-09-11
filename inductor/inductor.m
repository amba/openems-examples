#!/usr/bin/env -S octave --persist
                %
% Tutorials / MSL_NotchFilter
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Microstrip_Notch_Filter
%
% Tested with
%  - Matlab 2011a / Octave 4.0
%  - openEMS v0.0.33
%
% (C) 2011-2015 Thorsten Liebig <thorsten.liebig@gmx.de>

analysis_only = 0;

pkg load openems;
pkg load csxcad;


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
wire_radius = 100
num_loops = 10
loop_radius = 3e3
lead_height = 3 * loop_radius
coil_res = 100



%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD();
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = {'PML_8' 'PML_8' 'MUR' 'MUR' 'PEC' 'MUR'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
resolution = wire_radius/3;

mesh.x = SmoothMeshLines( [-MSL_length  MSL_length], resolution);
mesh.y = SmoothMeshLines( [-15*MSL_width 15*MSL_width], resolution);
mesh.z = SmoothMeshLines( [linspace(0,substrate_thickness,5), 2*lead_height], resolution );
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
portstop  = [ -separation/2+4*wire_radius,  MSL_width/2, 0];
[CSX,port{1}] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'FeedShift', 10*resolution, 'MeasPlaneShift',  MSL_length/3);

portstart = [mesh.x(end), -MSL_width/2, substrate_thickness];
portstop  = [separation/2-4*wire_radius         ,  MSL_width/2, 0];
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
count = count+1
p(1,count) = separation/2;
p(2,count) = 0;
p(3,count) = substrate_thickness;
CSX = AddWire(CSX, 'PEC', 0, p, wire_radius );


%% write/show/run the openEMS compatible xml-file
Sim_Path = 'tmp_notch_filter';
Sim_CSX = 'msl.xml';


if (analysis_only == 0)
  [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

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
