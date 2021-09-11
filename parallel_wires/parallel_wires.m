#!/usr/bin/env -S octave --persist

pkg load openems;
pkg load csxcad;


%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-6; % specify everything in um
MSL_length = 20e3;
MSL_width = 600;
substrate_thickness = 254;
substrate_epr = 3.66;

f_max = 7e9;

%% metal wires
wire_separation = 1000
wire_radius = 100
wire_height = 10e3



%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD();
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
BC   = {'PML_8' 'PML_8' 'MUR' 'MUR' 'PEC' 'MUR'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
resolution = c0/(f_max*sqrt(substrate_epr))/unit /50; % resolution of lambda/50
mesh.x = SmoothMeshLines( [0, 1.5*wire_separation], resolution/4, 1.5 ,0 );
mesh.x = SmoothMeshLines( [-MSL_length -mesh.x mesh.x MSL_length], resolution, 1.5 ,0 );
mesh.y = SmoothMeshLines( [0, MSL_width/2+[-resolution/3 +resolution/3*2]/4], resolution/4 , 1.5 ,0);
mesh.y = SmoothMeshLines( [-15*MSL_width, -mesh.y, mesh.y,  15*MSL_width], resolution/4);
mesh.z = SmoothMeshLines([0, substrate_thickness, 1.5*wire_height], resolution/4);

%% substrate
CSX = AddMaterial( CSX, 'RO4350B' );
CSX = SetMaterialProperty( CSX, 'RO4350B', 'Epsilon', substrate_epr );
start = [mesh.x(1),   mesh.y(1),   0];
stop  = [mesh.x(end), mesh.y(end), substrate_thickness];
CSX = AddBox( CSX, 'RO4350B', 0, start, stop );

CSX = AddMetal( CSX, 'PEC' );

%% wires
p(1,1) = -wire_separation/2;
p(2,1) = 0;
p(3,1) = substrate_thickness;

p(1,2) = -wire_separation/2;
p(2,2) = 0;
p(3,2) = substrate_thickness + wire_height;

CSX = AddWire(CSX, 'PEC', 0, p, wire_radius );

p(1,1) = wire_separation/2;
p(2,1) = 0;
p(3,1) = substrate_thickness;

p(1,2) = wire_separation/2;
p(2,2) = 0;
p(3,2) = substrate_thickness + wire_height;

CSX = AddWire(CSX, 'PEC', 0, p, wire_radius );

CSX = DefineRectGrid( CSX, unit, mesh );


%% MSL port

portstart = [ mesh.x(1), -MSL_width/2, substrate_thickness];
portstop  = [ -(wire_separation/2 - 2*wire_radius),  MSL_width/2, 0];
[CSX,port{1}] = AddMSLPort( CSX, 999, 1, 'PEC', portstart, portstop, 0, [0 0 -1], 'ExcitePort', true, 'FeedShift', 10*resolution, 'MeasPlaneShift',  MSL_length/3);

portstart = [mesh.x(end), -MSL_width/2, substrate_thickness];
portstop  = [(wire_separation/2 - 2*wire_radius),  MSL_width/2, 0];
[CSX,port{2}] = AddMSLPort( CSX, 999, 2, 'PEC', portstart, portstop, 0, [0 0 -1], 'MeasPlaneShift',  MSL_length/3 );




 
%% write/show/run the openEMS compatible xml-file
Sim_Path = 'tmp';
Sim_CSX = 'msl.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
RunOpenEMS( Sim_Path, Sim_CSX);

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
## ylim([-40 2]);
