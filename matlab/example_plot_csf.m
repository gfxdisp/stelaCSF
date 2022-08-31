% Plot CSF functions

csf_model = CSF_stelaCSF();

% 2D plot - as the function of temporal frequency
figure(1);
clf;

t_freq = linspace( 0, 60 )'; %Hz, must be a column vector
csf_pars = struct( 's_frequency', 4, 't_frequency', t_freq, 'orientation', 0, 'luminance', 100, 'area', 1, 'eccentricity', 0 );          
S = csf_model.sensitivity( csf_pars );        

plot( t_freq, S );
set( gca, 'YScale', 'log' );
xlabel( 'Temporal frequency [Hz]' );
ylabel( 'Sensitivity' );

% 3D plot - as the function of spatial and temporal frequency
figure(2)
clf;

t_freq = logspace( 0.1, log10(64), 30 );
s_freq = logspace( log10(0.5), log10(64), 30 );

[ss, tt] = meshgrid( s_freq, t_freq );

csf_pars = struct( 's_frequency', ss(:), 't_frequency', tt(:), 'orientation', 0, 'luminance', 100, 'area', 1, 'eccentricity', 0 );          
S = csf_model.sensitivity( csf_pars );        
S = reshape( S, size(ss) );

surf( s_freq, t_freq, S, 'FaceColor', 'interp', 'FaceLighting', 'phong' );
set( gca, 'XScale', 'log' );
set( gca, 'YScale', 'log' );
set( gca, 'ZScale', 'log' );
zlim( [1 1000] );
xlabel( 'Spatial frequency [cpd]')
ylabel( 'Temporal frequency [Hz]')
zlabel( 'Sensitivity')
title( 'stelaCSF');
