% MATLAB file must be saved as "Askisi_31_app_exported.m" for it to work properly.

classdef Askisi_31_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        GridLayout                 matlab.ui.container.GridLayout
        StartLabel                 matlab.ui.control.Label
        Peclet                     matlab.ui.control.Label
        IdealPecletLabel           matlab.ui.control.Label
        StepSpinner                matlab.ui.control.Spinner
        StepSpinnerLabel           matlab.ui.control.Label
        MaxPeSpinner               matlab.ui.control.Spinner
        MaxPeSpinnerLabel          matlab.ui.control.Label
        MinPeSpinner               matlab.ui.control.Spinner
        MinPeSpinnerLabel          matlab.ui.control.Label
        ConcentrationProfileLabel  matlab.ui.control.Label
        BetaSpinner                matlab.ui.control.Spinner
        Label                      matlab.ui.control.Label
        DaSpinner                  matlab.ui.control.Spinner
        DaSpinnerLabel             matlab.ui.control.Label
        BiSpinner                  matlab.ui.control.Spinner
        BiSpinnerLabel             matlab.ui.control.Label
        PeSpinner                  matlab.ui.control.Spinner
        PeSpinnerLabel             matlab.ui.control.Label
        NeSpinner                  matlab.ui.control.Spinner
        NeSpinnerLabel             matlab.ui.control.Label
        TemperatureProfileLabel    matlab.ui.control.Label
        ReactorOutAxis             matlab.ui.control.UIAxes
        ConcProfile                matlab.ui.control.UIAxes
        TempProfile                matlab.ui.control.UIAxes
        ConcGraph                  matlab.ui.control.UIAxes
        TempGraph                  matlab.ui.control.UIAxes
    end


    methods (Access = private)
        
        function[T, z] = SolveTemp(~, Ne, Pe, Bi)

            N = Ne + 1;
        
            A = zeros(N);
            b = zeros(N, 1);
            Z = zeros(Ne, 2);
        
            dz = 1. / Ne;
        
            Ml = dz / 6. * [
                2, 1 ;
                1, 2
            ];
            
            Ll = 1. / dz * [
                1, -1 ;
                -1, 1
            ];
            
            Cl = 0.5 * [
                -1, 1 ;
                -1, 1
            ];
        
            for ie = 1 : Ne
                Z(ie, 1) = ie;
                Z(ie, 2) = ie + 1;
            end
            
            for ie = 1 : Ne
                for i = 1 : 2
                    il = Z(ie, i);
            
                    for j = 1 : 2
                        jl = Z(ie, j);
            
                        A(il, jl) = A(il, jl) + Cl(i, j) + Ll(i, j) / Pe + Ml(i, j) * Bi / Pe;
                    end
                end
            end
            
            A(1, :) = 0;
            A(1, 1) = 1;
            b(1) = 1;
            
            b(N) = b(N) + 0;
            
            T = A \ b;
            z = 0 : dz : 1;
        end

        function [C, z] = SolveConc(~, Ne, Pe, Da, beta, T)
        
            N = Ne + 1;
        
            A = zeros(N);
            b = zeros(N, 1);
            Z = zeros(Ne, 2);
        
            dz = 1. / Ne;
        
            Ml = dz / 6. * [
                2, 1 ;
                1, 2
            ];
            
            Ll = 1. / dz * [
                1, -1 ;
                -1, 1
            ];
            
            Cl = 0.5 * [
                -1, 1 ;
                -1, 1
            ];
        
            for ie = 1 : Ne
                Z(ie, 1) = ie;
                Z(ie, 2) = ie + 1;
            end
            
            for ie = 1 : Ne
            
                T_e = (T(Z(ie, 1)) + T(Z(ie, 2))) / 2.;
        
                for i = 1 : 2
                    il = Z(ie, i);
            
                    for j = 1 : 2
                        jl = Z(ie, j);
            
                        A(il, jl) = A(il, jl) + Cl(i, j) + Ll(i, j) / Pe + exp(-beta * T_e) * Ml(i, j) * Da / Pe;
                    end
                end
            end
            
            A(1, :) = 0;
            A(1, 1) = 1;
            b(1) = 1;
            
            b(N) = b(N) + 0;
            
            C = A \ b;
            z = 0 : dz : 1;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: BetaSpinner, BiSpinner, DaSpinner, 
        % ...and 5 other components
        function Update(app, event)
            % Hide start label
            app.StartLabel.Visible = "off";
            
            % Get values
            Ne = app.NeSpinner.Value;
            Pe = app.PeSpinner.Value;
            Bi = app.BiSpinner.Value;
            Da = app.DaSpinner.Value;
            beta = app.BetaSpinner.Value;
            
            % Get profiles
            [T, z] = SolveTemp(app, Ne, Pe, Bi);
            c = SolveConc(app, Ne, Pe, Da, beta, T);

            % Plot temperature graph
            plot(app.TempGraph, z, T)
            xlim(app.TempGraph, [0, 1])
            ylim(app.TempGraph, [0, 1])

            % Plot concentration graph
            plot(app.ConcGraph, z, c)
            xlim(app.ConcGraph, [0, 1])
            ylim(app.ConcGraph, [0, 1])

            % 2D plot of temperature profile
            imagesc(app.TempProfile, z, [0, 1], T')
            clim(app.TempProfile, [0, 1])
            colorbar(app.TempProfile)

            % 2D plot of concentration profile
            imagesc(app.ConcProfile, z, [0, 1], c')
            clim(app.ConcProfile, [0, 1])
            colorbar(app.ConcProfile)

            % Reactor output
            PeMin = app.MinPeSpinner.Value;
            PeMax = app.MaxPeSpinner.Value;
            Step = app.StepSpinner.Value;

            if (PeMin >= PeMax)
                app.MaxPeSpinner.Value = PeMin + 20;
                PeMax = app.MaxPeSpinner.Value;
            end

            PeList = PeMin:Step:PeMax;

            Y = zeros(length(PeList), 1);
            for i = 1 : length(PeList)
                Temp = SolveTemp(app, Ne, PeList(i), Bi);
                Conc = SolveConc(app, Ne, PeList(i), Da, beta, Temp);

                c_exit = Conc(end);
                Y(i) = PeList(i) * (1 - c_exit);
            end

            PeIndex = Y == max(Y);
            app.Peclet.Text = sprintf("%d", PeList(PeIndex));

            plot(app.ReactorOutAxis, PeList, Y)
            xlim(app.ReactorOutAxis, [PeMin, PeMax])
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1111 708];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create TempGraph
            app.TempGraph = uiaxes(app.GridLayout);
            xlabel(app.TempGraph, 'z')
            ylabel(app.TempGraph, 'T')
            zlabel(app.TempGraph, 'Z')
            app.TempGraph.Box = 'on';
            app.TempGraph.XGrid = 'on';
            app.TempGraph.YGrid = 'on';
            app.TempGraph.Layout.Row = [2 6];
            app.TempGraph.Layout.Column = [2 6];

            % Create ConcGraph
            app.ConcGraph = uiaxes(app.GridLayout);
            xlabel(app.ConcGraph, 'z')
            ylabel(app.ConcGraph, 'c')
            zlabel(app.ConcGraph, 'Z')
            app.ConcGraph.Box = 'on';
            app.ConcGraph.XGrid = 'on';
            app.ConcGraph.YGrid = 'on';
            app.ConcGraph.Layout.Row = [2 6];
            app.ConcGraph.Layout.Column = [12 16];

            % Create TempProfile
            app.TempProfile = uiaxes(app.GridLayout);
            zlabel(app.TempProfile, 'Z')
            app.TempProfile.YTick = [];
            app.TempProfile.Box = 'on';
            app.TempProfile.Layout.Row = [7 8];
            app.TempProfile.Layout.Column = [2 6];

            % Create ConcProfile
            app.ConcProfile = uiaxes(app.GridLayout);
            zlabel(app.ConcProfile, 'Z')
            app.ConcProfile.YTick = [];
            app.ConcProfile.Box = 'on';
            app.ConcProfile.Layout.Row = [7 8];
            app.ConcProfile.Layout.Column = [12 16];

            % Create ReactorOutAxis
            app.ReactorOutAxis = uiaxes(app.GridLayout);
            title(app.ReactorOutAxis, 'Reactor Output')
            xlabel(app.ReactorOutAxis, 'Pe')
            ylabel(app.ReactorOutAxis, 'Pe(1-c)')
            zlabel(app.ReactorOutAxis, 'Z')
            app.ReactorOutAxis.FontName = 'Garamond';
            app.ReactorOutAxis.Box = 'on';
            app.ReactorOutAxis.XGrid = 'on';
            app.ReactorOutAxis.YGrid = 'on';
            app.ReactorOutAxis.FontSize = 20;
            app.ReactorOutAxis.Layout.Row = [8 12];
            app.ReactorOutAxis.Layout.Column = [7 11];

            % Create TemperatureProfileLabel
            app.TemperatureProfileLabel = uilabel(app.GridLayout);
            app.TemperatureProfileLabel.HorizontalAlignment = 'center';
            app.TemperatureProfileLabel.VerticalAlignment = 'bottom';
            app.TemperatureProfileLabel.FontName = 'Garamond';
            app.TemperatureProfileLabel.FontSize = 20;
            app.TemperatureProfileLabel.FontWeight = 'bold';
            app.TemperatureProfileLabel.Layout.Row = 1;
            app.TemperatureProfileLabel.Layout.Column = [3 5];
            app.TemperatureProfileLabel.Text = 'Temperature Profile';

            % Create NeSpinnerLabel
            app.NeSpinnerLabel = uilabel(app.GridLayout);
            app.NeSpinnerLabel.HorizontalAlignment = 'right';
            app.NeSpinnerLabel.FontName = 'Garamond';
            app.NeSpinnerLabel.FontSize = 20;
            app.NeSpinnerLabel.Layout.Row = 2;
            app.NeSpinnerLabel.Layout.Column = 8;
            app.NeSpinnerLabel.Text = 'Ne';

            % Create NeSpinner
            app.NeSpinner = uispinner(app.GridLayout);
            app.NeSpinner.Limits = [2 Inf];
            app.NeSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.NeSpinner.HorizontalAlignment = 'center';
            app.NeSpinner.FontName = 'Garamond';
            app.NeSpinner.FontSize = 20;
            app.NeSpinner.Layout.Row = 2;
            app.NeSpinner.Layout.Column = [9 10];
            app.NeSpinner.Value = 20;

            % Create PeSpinnerLabel
            app.PeSpinnerLabel = uilabel(app.GridLayout);
            app.PeSpinnerLabel.HorizontalAlignment = 'right';
            app.PeSpinnerLabel.FontName = 'Garamond';
            app.PeSpinnerLabel.FontSize = 20;
            app.PeSpinnerLabel.Layout.Row = 3;
            app.PeSpinnerLabel.Layout.Column = 8;
            app.PeSpinnerLabel.Text = 'Pe';

            % Create PeSpinner
            app.PeSpinner = uispinner(app.GridLayout);
            app.PeSpinner.Limits = [1 Inf];
            app.PeSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.PeSpinner.HorizontalAlignment = 'center';
            app.PeSpinner.FontName = 'Garamond';
            app.PeSpinner.FontSize = 20;
            app.PeSpinner.Layout.Row = 3;
            app.PeSpinner.Layout.Column = [9 10];
            app.PeSpinner.Value = 20;

            % Create BiSpinnerLabel
            app.BiSpinnerLabel = uilabel(app.GridLayout);
            app.BiSpinnerLabel.HorizontalAlignment = 'right';
            app.BiSpinnerLabel.FontName = 'Garamond';
            app.BiSpinnerLabel.FontSize = 20;
            app.BiSpinnerLabel.Layout.Row = 4;
            app.BiSpinnerLabel.Layout.Column = 8;
            app.BiSpinnerLabel.Text = 'Bi';

            % Create BiSpinner
            app.BiSpinner = uispinner(app.GridLayout);
            app.BiSpinner.Limits = [1 Inf];
            app.BiSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.BiSpinner.HorizontalAlignment = 'center';
            app.BiSpinner.FontName = 'Garamond';
            app.BiSpinner.FontSize = 20;
            app.BiSpinner.Layout.Row = 4;
            app.BiSpinner.Layout.Column = [9 10];
            app.BiSpinner.Value = 50;

            % Create DaSpinnerLabel
            app.DaSpinnerLabel = uilabel(app.GridLayout);
            app.DaSpinnerLabel.HorizontalAlignment = 'right';
            app.DaSpinnerLabel.FontName = 'Garamond';
            app.DaSpinnerLabel.FontSize = 20;
            app.DaSpinnerLabel.Layout.Row = 5;
            app.DaSpinnerLabel.Layout.Column = 8;
            app.DaSpinnerLabel.Text = 'Da';

            % Create DaSpinner
            app.DaSpinner = uispinner(app.GridLayout);
            app.DaSpinner.Limits = [1 Inf];
            app.DaSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.DaSpinner.HorizontalAlignment = 'center';
            app.DaSpinner.FontName = 'Garamond';
            app.DaSpinner.FontSize = 20;
            app.DaSpinner.Layout.Row = 5;
            app.DaSpinner.Layout.Column = [9 10];
            app.DaSpinner.Value = 1000;

            % Create Label
            app.Label = uilabel(app.GridLayout);
            app.Label.HorizontalAlignment = 'right';
            app.Label.FontName = 'Garamond';
            app.Label.FontSize = 20;
            app.Label.Layout.Row = 6;
            app.Label.Layout.Column = 8;
            app.Label.Text = 'Beta';

            % Create BetaSpinner
            app.BetaSpinner = uispinner(app.GridLayout);
            app.BetaSpinner.Limits = [0 Inf];
            app.BetaSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.BetaSpinner.HorizontalAlignment = 'center';
            app.BetaSpinner.FontName = 'Garamond';
            app.BetaSpinner.FontSize = 20;
            app.BetaSpinner.Layout.Row = 6;
            app.BetaSpinner.Layout.Column = [9 10];
            app.BetaSpinner.Value = 10;

            % Create ConcentrationProfileLabel
            app.ConcentrationProfileLabel = uilabel(app.GridLayout);
            app.ConcentrationProfileLabel.HorizontalAlignment = 'center';
            app.ConcentrationProfileLabel.VerticalAlignment = 'bottom';
            app.ConcentrationProfileLabel.FontName = 'Garamond';
            app.ConcentrationProfileLabel.FontSize = 20;
            app.ConcentrationProfileLabel.FontWeight = 'bold';
            app.ConcentrationProfileLabel.Layout.Row = 1;
            app.ConcentrationProfileLabel.Layout.Column = [13 15];
            app.ConcentrationProfileLabel.Text = 'Concentration Profile';

            % Create MinPeSpinnerLabel
            app.MinPeSpinnerLabel = uilabel(app.GridLayout);
            app.MinPeSpinnerLabel.HorizontalAlignment = 'right';
            app.MinPeSpinnerLabel.FontName = 'Garamond';
            app.MinPeSpinnerLabel.FontSize = 20;
            app.MinPeSpinnerLabel.Layout.Row = 9;
            app.MinPeSpinnerLabel.Layout.Column = [11 12];
            app.MinPeSpinnerLabel.Text = 'Min. Pe';

            % Create MinPeSpinner
            app.MinPeSpinner = uispinner(app.GridLayout);
            app.MinPeSpinner.Limits = [0 Inf];
            app.MinPeSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.MinPeSpinner.HorizontalAlignment = 'center';
            app.MinPeSpinner.FontName = 'Garamond';
            app.MinPeSpinner.FontSize = 20;
            app.MinPeSpinner.Layout.Row = 9;
            app.MinPeSpinner.Layout.Column = [13 14];
            app.MinPeSpinner.Value = 10;

            % Create MaxPeSpinnerLabel
            app.MaxPeSpinnerLabel = uilabel(app.GridLayout);
            app.MaxPeSpinnerLabel.HorizontalAlignment = 'right';
            app.MaxPeSpinnerLabel.FontName = 'Garamond';
            app.MaxPeSpinnerLabel.FontSize = 20;
            app.MaxPeSpinnerLabel.Layout.Row = 10;
            app.MaxPeSpinnerLabel.Layout.Column = [11 12];
            app.MaxPeSpinnerLabel.Text = 'Max. Pe';

            % Create MaxPeSpinner
            app.MaxPeSpinner = uispinner(app.GridLayout);
            app.MaxPeSpinner.Limits = [0 Inf];
            app.MaxPeSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.MaxPeSpinner.HorizontalAlignment = 'center';
            app.MaxPeSpinner.FontName = 'Garamond';
            app.MaxPeSpinner.FontSize = 20;
            app.MaxPeSpinner.Layout.Row = 10;
            app.MaxPeSpinner.Layout.Column = [13 14];
            app.MaxPeSpinner.Value = 200;

            % Create StepSpinnerLabel
            app.StepSpinnerLabel = uilabel(app.GridLayout);
            app.StepSpinnerLabel.HorizontalAlignment = 'right';
            app.StepSpinnerLabel.FontName = 'Garamond';
            app.StepSpinnerLabel.FontSize = 20;
            app.StepSpinnerLabel.Layout.Row = 11;
            app.StepSpinnerLabel.Layout.Column = 12;
            app.StepSpinnerLabel.Text = 'Step';

            % Create StepSpinner
            app.StepSpinner = uispinner(app.GridLayout);
            app.StepSpinner.Limits = [1 Inf];
            app.StepSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.StepSpinner.HorizontalAlignment = 'center';
            app.StepSpinner.FontName = 'Garamond';
            app.StepSpinner.FontSize = 20;
            app.StepSpinner.Layout.Row = 11;
            app.StepSpinner.Layout.Column = [13 14];
            app.StepSpinner.Value = 10;

            % Create IdealPecletLabel
            app.IdealPecletLabel = uilabel(app.GridLayout);
            app.IdealPecletLabel.HorizontalAlignment = 'center';
            app.IdealPecletLabel.VerticalAlignment = 'bottom';
            app.IdealPecletLabel.FontName = 'Garamond';
            app.IdealPecletLabel.FontSize = 20;
            app.IdealPecletLabel.FontWeight = 'bold';
            app.IdealPecletLabel.Layout.Row = 9;
            app.IdealPecletLabel.Layout.Column = [4 6];
            app.IdealPecletLabel.Text = 'Ideal Peclet';

            % Create Peclet
            app.Peclet = uilabel(app.GridLayout);
            app.Peclet.HorizontalAlignment = 'center';
            app.Peclet.VerticalAlignment = 'top';
            app.Peclet.FontName = 'Garamond';
            app.Peclet.FontSize = 20;
            app.Peclet.Layout.Row = 10;
            app.Peclet.Layout.Column = 5;
            app.Peclet.Text = '0';

            % Create StartLabel
            app.StartLabel = uilabel(app.GridLayout);
            app.StartLabel.HorizontalAlignment = 'center';
            app.StartLabel.FontName = 'Garamond';
            app.StartLabel.FontSize = 20;
            app.StartLabel.FontWeight = 'bold';
            app.StartLabel.Layout.Row = 1;
            app.StartLabel.Layout.Column = [7 12];
            app.StartLabel.Text = 'Modify any value to start plotting';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Askisi_31_app_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
