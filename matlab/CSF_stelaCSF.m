classdef CSF_stelaCSF < CSF_base
    % The class predicts contrast sensitivity as the functioon of
    % Spatio-Temporal frequency, Eccentricity, Luminance and Area (thus
    % stela).
    %
    % The details regarding the model, its derivation and calibration can
    % be found in the paper:
    %
    % Mantiuk, Rafał K, Maliha Ashraf, and Alexandre Chapiro.
    % “StelaCSF - A Unified Model of Contrast Sensitivity as the Function of Spatio-Temporal Frequency , Eccentricity , Luminance and Area.”
    % ACM Transactions on Graphics 41, no. 4 (2022): 145.
    % https://doi.org/10.1145/3528223.3530115.

    properties( Constant )
        Y_min = 0.001;  % The minimum luminance
        Y_max = 10000;  % The maximum luminance
        rho_min = 2^-4  % The minimum spatial frequency
        rho_max = 64;   % The maximum spatial frequency
        ecc_max = 120;  % The maximum eccentricity
    end

    properties
        use_gpu = true;
        ps_beta = 1;
    end

    methods

        function obj = CSF_stelaCSF( )
            obj.par = obj.get_default_par();
        end

        function name = short_name( obj )
            % A short name that could be used as a part of a file name
            name = 'stela-csf';
        end

        function name = full_name( obj )
            name = 'stelaCSF';
        end

        % Return contrast sensitivity for a given set of parameters. 
        % The sensitivity is assumed to be either the inverse of luminance
        % cotrast (L/\Delta L) or the inverse of cone contrast.
        %
        % This is a general interface that takes as input a
        % structure with the parameters values. This allows to add/remove
        % parameters without changing the interface, more flexible use of
        % parameters (e.g. either LMS or luminance of the background) and
        % should be less error prone. 
        % 
        % pars is the structure with the field names that closely match
        % the column names in the data files:
        %
        % luminance - luminance in cd/m^2. D65 background is assumed. 
        % lms_bkg - specify background colour and luminance 
        % s_frequency - spatial frequency in cpd
        % t_frequency - temporal frequency in Hz (default is 0)
        % orientation - orientation in deg (default is 0)
        % lms_delta - modulation direction in the LMS colour space (default
        %             D65 luminance modulation)
        % area - area of the stimulus in deg^2
        % ge_sigma - radius of the gaussian envelope in deg 
        % eccentricity - eccentricity in deg (default 0)       
        % vis_field - orientation in the visual field. See README.md
        %             (default 0)
        %
        % 'lms_bkg' and 'lms_delta' should be either [1 3] or [N 3] matrix. 
        % Any other field should be a column vector of size N or a scalar. 
        % For best performance, pass vectors with the the required number
        % of parameters. Do not call the sensitibity() function repetitively.
        % 
        %
        % You must specify 'luminance' or 'lms_bkg' but not both. 
        % You must specify 'area' or 'ge_sigma' but not both. 
        %
        % Example: 
        % csf_model = CSF_stelaCSF();
        % csf_pars = struct( 's_frequency', 4, 't_frequency', 1, 'orientation', 0, 'lms_bkg', [0.7443 0.3054 0.0157], 'area', 1, 'eccentricity', 0 );          
        % S = csf_model.sensitivity( csf_pars );        
        function S = sensitivity( obj, csf_pars )

            csf_pars = obj.test_complete_params(csf_pars, { 'luminance', 'ge_sigma' } );

            ecc = csf_pars.eccentricity;
            sigma = csf_pars.ge_sigma;
            rho = csf_pars.s_frequency;
            omega = csf_pars.t_frequency;
            lum = csf_pars.luminance;

            [R_sust, R_trans] = get_sust_trans_resp(obj, omega);

            A = pi*(sigma).^2; % Stimulus area

            S_sust = obj.csf_achrom( rho, A, lum, ecc, obj.par.ach_sust );
            S_trans = obj.csf_achrom( rho, A, lum, ecc, obj.par.ach_trans );

%             if obj.ps_beta ~= 1
%                 beta = obj.ps_beta;
%                 S = ( (R_sust.*S_sust).^beta + (R_trans.*S_trans).^beta).^(1/beta);
%             else
%                 S = R_sust.*S_sust + R_trans.*S_trans;
%             end

            S_aux = 0; %obj.aux_sensitivity( csf_pars );
            pm_ratio=1;
            if obj.ps_beta ~= 1 
                beta = obj.ps_beta;
                S = ( (R_sust.*S_sust .* sqrt(pm_ratio)).^beta + (R_trans.*S_trans .* sqrt(1./pm_ratio)).^beta + S_aux.^beta).^(1/beta);
            else
                S = R_sust.*S_sust .* sqrt(pm_ratio) + R_trans.*S_trans .* sqrt(1./pm_ratio) + S_aux;
            end
            

            % The drop of sensitivity with the eccentricity (the window of
            % visibiliy model + extension)
            alpha = min(1, abs(csf_pars.vis_field-180)/90 );
            ecc_drop = alpha .* obj.par.ecc_drop + (1-alpha) .* obj.par.ecc_drop_nasal;
            ecc_drop_f = alpha .* obj.par.ecc_drop_f + (1-alpha) .* obj.par.ecc_drop_f_nasal;
            a = ecc_drop + rho.*ecc_drop_f;
            S = S .* 10.^(-a.*ecc);

        end

        % Get the sustained and transient temporal response functions
        % omega - temporal frequency in Hz
        function [R_sust, R_trans] = get_sust_trans_resp(obj, omega)
            sigma_sust = obj.par.sigma_sust;
            beta_sust = 1.3314;
            omega_0 = 5;

            beta_trans = 0.1898;
            sigma_trans = obj.par.sigma_trans;

            R_sust = exp( -omega.^beta_sust / (sigma_sust) );
            R_trans = exp( -abs(omega.^beta_trans-omega_0^beta_trans).^2 / (sigma_trans) );
        end

        % Achromatic CSF model
        function S = csf_achrom( obj, freq, area, lum, ecc, ach_pars )
            % Internal. Do not call from outside the object.
            % A nested CSF as a function of luminance

            N = max( [length(freq) length(area) length(lum)] );

            assert( length(freq)==1 || all( size(freq)==[N 1] ) );
            assert( length(area)==1 || all( size(area)==[N 1] ) );
            assert( length(lum)==1 || all( size(lum)==[N 1] ) );

            S_max = obj.get_lum_dep( ach_pars.S_max, lum );
            f_max = obj.get_lum_dep( ach_pars.f_max, lum );
            %gamma = obj.get_lum_dep( ach_pars.gamma, lum );

            bw = ach_pars.bw;
            %bw = max( 0.01, ach_pars.bw - ach_pars.ecc_bw_drop*ecc);
            a = ach_pars.a;

            % Truncated log-parabola
            S_LP = 10.^( -(log10(freq) - log10(f_max)).^2./(2.^bw) );
            ss = (freq<f_max) & (S_LP < (1-a));
            S_LP(ss) = 1-a;

            S_peak = S_max .* S_LP;


            % The stimulus size model from the paper:
            %
            % Rovamo, J., Luntinen, O., & N�s�nen, R. (1993).
            % Modelling the dependence of contrast sensitivity on grating area and spatial frequency.
            % Vision Research, 33(18), 2773�2788.
            %
            % Equation on the page 2784, one after (25)

            if isfield( ach_pars, 'f0' )
                f0 = ach_pars.f0;
            else
                f0 = 0.65;
            end
            if isfield( ach_pars, 'A0' )
                A0 = ach_pars.A0;
            else
                A0 = 270;
            end

            Ac = A0./(1+(freq/f0).^2);

            S = S_peak .* sqrt( Ac ./ (1+Ac./area)).*(freq.^1);
        end

        function pd = get_plot_description( obj )
            pd = struct();
            pp = 1;
            pd(pp).title = 'Sustained and transient response';
            pd(pp).id = 'sust_trans';
            pp = pp+1;
            pd(pp).title = 'Peak sensitivity';
            pd(pp).id = 'peak_s';
            pp = pp+1;
        end

        function plot_mechanism( obj, plt_id )
            switch( plt_id )
                case 'sust_trans' % sust-trans-response
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    omega = linspace( 0, 100 );
                    [R_sust, R_trans] = obj.get_sust_trans_resp(omega);
                    hh(1) = plot( omega, R_sust, 'DisplayName', 'Sustained');
                    hold on
                    hh(2) = plot( omega, R_trans, 'DisplayName', 'Transient');
                    hold off
                    xlabel( 'Temp. freq. [Hz]' );
                    ylabel( 'Response' );
                    legend( hh, 'Location', 'Best' );
                    grid on;

                case { 'peak_s' }

                    f = logspace( -2, log10(5), 1024 );
                    L = logspace( -2, 4 );
                    [LL, ff] = meshgrid( L, f );
                    OMEGAs = [0 5 16];
                    COLORs = lines(length(OMEGAs));

                    for pp=1:length(OMEGAs)
                        csfpar.luminance = LL(:);
                        csfpar.s_frequency = ff(:);
                        csfpar.t_frequency = OMEGAs(pp);
                        %csfpar.area = pi*(0.5./ff(:)).^2;
                        csfpar.area = pi*(1.5).^2;
                        csfpar.eccentricity = 0;

                        S = obj.sensitivity( csfpar );

                        S = reshape( S, size(ff) );

                        S_max = max(S);

                        hh(pp) = plot( L, S_max, 'Color', COLORs(pp,:), 'DisplayName', sprintf('%g Hz', OMEGAs(pp)) );
                        hold on
                    end

                    L_dvr = logspace( -1, 0 );
                    hh(pp+1) = plot( L_dvr, sqrt(L_dvr)*50, '--k', 'DisplayName', 'DeVries-Rose law' );

                    set_axis_tick_label( 'x', 'luminance', L );
                    set_axis_tick_label( 'y', 'sensitivity', [1 1000]);

                    legend( hh, 'Location', 'best' );
                    ylabel( 'Peak sensitivity' );

                    ylim( [1 1000] );
                    grid on


                case { 'sust_peak_s', 'trans_peak_s' }
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    L = logspace( -2, 4 );
                    if strcmp( plt_id, 'sust_peak_s' )
                        S_max = obj.par.ach_sust.S_max;
                    else
                        S_max = obj.par.ach_trans.S_max;
                    end
                    plot( L,  obj.get_lum_dep( S_max, L ) );
                    hold on
                    L_dvr = logspace( -1, 1 );
                    plot( L_dvr, sqrt(L_dvr)*100, '--k' );

                    set_axis_tick_label( 'x', 'luminance', L );
                    set_axis_tick_label( 'y', 'sensitivity', [1 100000] );
                    grid on;
                case 'sust_peak_f'
                    clf;
                    html_change_figure_print_size( gcf, 10, 10 );
                    L = logspace( -2, 4 );
                    %                     if plt_id == 5
                    f_max = obj.par.ach_sust.f_max;
                    %                     else
                    %                         f_max = obj.par.ach_trans.f_max;
                    %                     end
                    plot( L,  obj.get_lum_dep( f_max, L ) );
                    set_axis_tick_label( 'x', 'luminance', L );
                    set_axis_tick_label( 'y', 'frequency', [0.01 60] );
                    grid on;
                otherwise
                    error( 'Wrong plt_id' );
            end
        end

    end

    methods( Static )

        function p = get_default_par()

            p = CSF_base.get_dataset_par();

            % fitting based on the run: stela-csf_all_2022-04-16_16-48
            p.ach_sust.S_max = [ 68.9501 59.5023 0.164274 7.54866e-07 7.77268e+09 ];
            p.ach_sust.f_max = [ 1.62144 36.6565 0.255823 ];
            p.ach_sust.bw = 0.000219263;
            p.ach_sust.a = 0.103686;
            p.ach_trans.S_max = [ 0.500846 57.3469 ];
            p.ach_trans.f_max = 0.0267489;
            p.ach_trans.bw = 1.75147;
            p.ach_trans.a = 0.000273289;
            p.sigma_trans = 0.12314;
            p.sigma_sust = 5.79336;
            p.ecc_drop = 0.0296662;
            p.ecc_drop_nasal = 0.0113638;
            p.ecc_drop_f = 0.0190062;
            p.ecc_drop_f_nasal = 0.0193858;

        end


    end

end