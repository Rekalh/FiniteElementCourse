% MATLAB file must be named "Askisi_32_app_exported" for it to work properly.

classdef Askisi_32_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        GridLayout           matlab.ui.container.GridLayout
        ForceNexNeyCheckBox  matlab.ui.control.CheckBox
        ModifyanyvaluetostartplottingLabel  matlab.ui.control.Label
        NeySpinner           matlab.ui.control.Spinner
        NeySpinnerLabel      matlab.ui.control.Label
        NexSpinner           matlab.ui.control.Spinner
        NexSpinnerLabel      matlab.ui.control.Label
        UIAxes               matlab.ui.control.UIAxes
    end

    
    methods (Access = private)
        
        function [T_ij, x, y] = Solve(~, Nex, Ney)
            Ne = Nex * Ney;

            Nx = Nex + 1;
            Ny = Ney + 1;
            N = Nx * Ny;
            
            A = zeros(N);
            b = zeros(N, 1);
            s = zeros(N, 1);
            Z = zeros(Ne, 4);
            
            dx = 2. / Nex;
            dy = 2. / Ney;
            
            for iex = 1 : Nex
                for iey = 1 : Ney
                    ie = iey + Ney * (iex - 1);
            
                    Z(ie, 1) = iey + Ny * (iex - 1);
                    Z(ie, 2) = Z(ie, 1) + 1;
                    Z(ie, 3) = Z(ie, 1) + Ny;
                    Z(ie, 4) = Z(ie, 3) + 1;
                end
            end
            
            for i = 1 : Nx
                for j = 1 : Ny
                    k = j + (i - 1) * Ny;
                    x = -1 + (i - 1) * dx;
                    y = -1 + (j - 1) * dy;
                    
                    s(k) = cos(pi * x / 2.) * cos(pi * y / 2.);
                end
            end
            
            Ml = [
                4., 2., 2., 1. ;
                2., 4., 1., 2. ;
                2., 1., 4., 2. ;
                1., 2., 2., 4.
            ];
            Ml = Ml * dx * dy / 36.;
            
            Llx = [
                2, 1, -2, -1 ;
                1, 2, -1, -2 ;
                -2, -1, 2, 1 ;
                -1, -2, 1, 2
            ];
            Llx = Llx * dy / dx / 6.;
            
            Lly = [
                2, -2, 1, -1 ;
                -2, 2, -1, 1 ;
                1, -1, 2, -2 ;
                -1, 1, -2, 2
            ];
            Lly = Lly * dx / dy / 6.;
            
            Ll = Llx + Lly;
            
            iflagdir = zeros(N, 1);
            
            for i = 1 : Nx
                iflagdir(1 + (i - 1) * Ny) = 1;
                iflagdir(Ny + (i - 1) * Ny) = 1;
            end
            
            for j = 1 : Ny
                iflagdir(j) = 1;
                iflagdir(j + (Nx - 1) * Ny) = 1;
            end
            
            for ie = 1 : Ne
                for i = 1 : 4
                    il = Z(ie, i);
            
                    for j = 1 : 4
                        jl = Z(ie, j);
            
                        A(il, jl) = A(il, jl) + Ll(i, j);
                        b(il) = b(il) + Ml(i, j) * s(jl);
                    end
                end
            end
            
            for i = 1 : N
                if (iflagdir(i) ~= 0)
                    A(i, :) = 0.;
                    A(i, i) = 1.;
                    b(i) = 1.;
                end
            end
            
            T = A \ b;
            
            T_ij = reshape(T, [Ny, Nx]);
            
            x = -1 : dx : 1;
            y = -1 : dy : 1;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: ForceNexNeyCheckBox, NexSpinner, 
        % ...and 1 other component
        function Update(app, event)
            Nex = app.NexSpinner.Value;
            Ney = app.NeySpinner.Value;
            ForceEqual = app.ForceNexNeyCheckBox.Value;

            if (ForceEqual)
                Nex = max(Nex, Ney);
                Ney = Nex;

                app.NexSpinner.Value = Nex;
                app.NeySpinner.Value = Nex;
            end

            app.ModifyanyvaluetostartplottingLabel.Text = "Solving...";
            drawnow
            
            [T_ij, x, y] = Solve(app, Nex, Ney);

            imagesc(app.UIAxes, x, y, T_ij)
            xlim(app.UIAxes, [-1, 1])
            ylim(app.UIAxes, [-1, 1])
            title(app.UIAxes, sprintf("Ne = %d", Nex * Ney))
            colorbar(app.UIAxes)

            app.ModifyanyvaluetostartplottingLabel.Text = "Solved";
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1000 695];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.XTick = [];
            app.UIAxes.XTickLabel = '';
            app.UIAxes.YTick = [];
            app.UIAxes.YTickLabel = '';
            app.UIAxes.Box = 'on';
            app.UIAxes.Layout.Row = [2 8];
            app.UIAxes.Layout.Column = [4 10];

            % Create NexSpinnerLabel
            app.NexSpinnerLabel = uilabel(app.GridLayout);
            app.NexSpinnerLabel.HorizontalAlignment = 'right';
            app.NexSpinnerLabel.FontName = 'Garamond';
            app.NexSpinnerLabel.FontSize = 20;
            app.NexSpinnerLabel.Layout.Row = 4;
            app.NexSpinnerLabel.Layout.Column = 1;
            app.NexSpinnerLabel.Text = 'Nex';

            % Create NexSpinner
            app.NexSpinner = uispinner(app.GridLayout);
            app.NexSpinner.Limits = [2 Inf];
            app.NexSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.NexSpinner.HorizontalAlignment = 'center';
            app.NexSpinner.FontName = 'Garamond';
            app.NexSpinner.FontSize = 20;
            app.NexSpinner.Layout.Row = 4;
            app.NexSpinner.Layout.Column = [2 3];
            app.NexSpinner.Value = 2;

            % Create NeySpinnerLabel
            app.NeySpinnerLabel = uilabel(app.GridLayout);
            app.NeySpinnerLabel.HorizontalAlignment = 'right';
            app.NeySpinnerLabel.FontName = 'Garamond';
            app.NeySpinnerLabel.FontSize = 20;
            app.NeySpinnerLabel.Layout.Row = 5;
            app.NeySpinnerLabel.Layout.Column = 1;
            app.NeySpinnerLabel.Text = 'Ney';

            % Create NeySpinner
            app.NeySpinner = uispinner(app.GridLayout);
            app.NeySpinner.Limits = [2 Inf];
            app.NeySpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.NeySpinner.HorizontalAlignment = 'center';
            app.NeySpinner.FontName = 'Garamond';
            app.NeySpinner.FontSize = 20;
            app.NeySpinner.Layout.Row = 5;
            app.NeySpinner.Layout.Column = [2 3];
            app.NeySpinner.Value = 2;

            % Create ModifyanyvaluetostartplottingLabel
            app.ModifyanyvaluetostartplottingLabel = uilabel(app.GridLayout);
            app.ModifyanyvaluetostartplottingLabel.HorizontalAlignment = 'center';
            app.ModifyanyvaluetostartplottingLabel.FontName = 'Garamond';
            app.ModifyanyvaluetostartplottingLabel.FontSize = 20;
            app.ModifyanyvaluetostartplottingLabel.Layout.Row = 9;
            app.ModifyanyvaluetostartplottingLabel.Layout.Column = [5 9];
            app.ModifyanyvaluetostartplottingLabel.Text = 'Modify any value to start plotting';

            % Create ForceNexNeyCheckBox
            app.ForceNexNeyCheckBox = uicheckbox(app.GridLayout);
            app.ForceNexNeyCheckBox.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.ForceNexNeyCheckBox.Text = 'Force Nex = Ney';
            app.ForceNexNeyCheckBox.FontName = 'Garamond';
            app.ForceNexNeyCheckBox.FontSize = 15;
            app.ForceNexNeyCheckBox.Layout.Row = 6;
            app.ForceNexNeyCheckBox.Layout.Column = [2 3];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Askisi_32_app_exported

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
