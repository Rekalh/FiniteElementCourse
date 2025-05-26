classdef Askisi_23_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        GridLayout                  matlab.ui.container.GridLayout
        InterpolationDropDown       matlab.ui.control.DropDown
        InterpolationDropDownLabel  matlab.ui.control.Label
        StopButton                  matlab.ui.control.Button
        StartButton                 matlab.ui.control.Button
        ShowGridCheckBox            matlab.ui.control.CheckBox
        ShowNumericalCheckBox       matlab.ui.control.CheckBox
        ShowAnalyticalCheckBox      matlab.ui.control.CheckBox
        tmaxEditField               matlab.ui.control.NumericEditField
        tmaxEditFieldLabel          matlab.ui.control.Label
        OmegaEditField              matlab.ui.control.NumericEditField
        OmegaEditFieldLabel         matlab.ui.control.Label
        LEditField                  matlab.ui.control.NumericEditField
        LEditFieldLabel             matlab.ui.control.Label
        TimeStepEditField           matlab.ui.control.NumericEditField
        TimeStepEditFieldLabel      matlab.ui.control.Label
        NSpinner                    matlab.ui.control.Spinner
        NSpinnerLabel               matlab.ui.control.Label
        UIAxes4                     matlab.ui.control.UIAxes
        UIAxes3                     matlab.ui.control.UIAxes
        UIAxes                      matlab.ui.control.UIAxes
        UIAxes2                     matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        s1
        s2

        AnimationTimer
        CurrentFrame = 1
        IsAnimating = false
        AnimationIterations = 0

        AnimationPlot
        AnimationPlotA
        Interpolation = "None";
    end
    
    methods (Access = private)
        
        function [x, tm, T, Ta] = FE_Askisi_23(~, N, dt, L, omega, tmax)

            Ne = N - 1;
            
            dx = L / N;
            
            Nt = round(tmax / dt - 1);
            
            k = sqrt(omega / 2);
            
            % Define matrices
            Ml = zeros(2);
            Ll = zeros(2);
            A = zeros(N);
            
            b = zeros(N, 1);
            Z = zeros(N, 2);
            tm = zeros(N, 1);
            T = zeros(N, Nt);
            Ta = zeros(N, Nt);
            x = zeros(N, 1);
            Tp = zeros(N, 1);
            
            % Fill connection matrix
            for ie = 1 : Ne
                Z(ie, 1) = ie;
                Z(ie, 2) = ie + 1;
            end
            
            % Fill x
            for i = 1 : N
                x(i) = (i - 1) * dx;
            end
            
            % Fill local matrices
            Ml(1, 1) = 2.;
            Ml(1, 2) = 1.;
            Ml(2, 1) = 1.;
            Ml(2, 2) = 2.;
            
            Ml = Ml * dx / 6;
            
            Ll(1, 1) = 1.;
            Ll(1, 2) = -1.;
            Ll(2, 1) = -1.;
            Ll(2, 2) = 1.;
            
            Ll = Ll / dx;
            
            % Initial condition
            T(1:N) = 0.;
            
            t = 0;
            
            c = 1;
            while (t < tmax)
                
                A(:, :) = 0;
                b(:) = 0;

                for i = 1 : N
                    Ta(i, c) = exp(-k * x(i)) * sin(omega * t - k * x(i));
                end
            
                for ie = 1 : Ne
                    for i = 1 : 2
                        il = Z(ie, i);
                        for j = 1 : 2
                            jl = Z(ie, j);
            
                            A(il, jl) = A(il, jl) + 2 * Ml(i, j) + Ll(i, j) * dt;
                            b(il) = b(il) + (2 * Ml(i, j) - Ll(i, j) * dt) * Tp(jl);
                        end
                    end
                end
            
                % Dirichlet boundary conditions
                A(1, 1:N) = 0;
                A(1, 1) = 1;
                b(1) = sin(omega * t);
            
                A(N, 1:N) = 0;
                A(N, N) = 1;
                b(N) = 0;
            
                T(1:N, c) = A \ b;
                Tp = T(1:N, c);
                tm(c) = t;
            
                t = t + dt;
                c = c + 1;
            end
        end
        
        function [] = updateAnimation(app)
            
            if (~app.IsAnimating)
                stop(app.AnimationTimer)
                return;
            end

            N = app.NSpinner.Value;
            dt = app.TimeStepEditField.Value;
            L = app.LEditField.Value;
            omega = app.OmegaEditField.Value;

            period = 2 * pi / omega;
            
            [x, t, T, Ta] = FE_Askisi_23(app, N, dt, L, omega, 15 * period);

            if (app.CurrentFrame > length(t))
                app.CurrentFrame = 1;
                app.AnimationIterations = app.AnimationIterations + 1;
            end

            Tq = T(:, app.CurrentFrame);
            xq = x;
            Taq = Ta(:, app.CurrentFrame);

            if (app.Interpolation ~= "None")
                xq = linspace(0, L, 100);
                Tq = interp1(x, T(:, app.CurrentFrame), xq, app.Interpolation)';
                Taq = interp1(x, Ta(:, app.CurrentFrame), xq, app.Interpolation);
            end

            imagesc(app.UIAxes4, xq, [0 1], Tq');
            clim(app.UIAxes4, [-1, 1]);
            xlim(app.UIAxes4, [0, L]);
            colorbar(app.UIAxes4);

            app.AnimationPlot.Visible = "on";
            app.AnimationPlot.YData = Tq;
            app.AnimationPlot.XData = xq;

            app.AnimationPlotA.Visible = "on";
            app.AnimationPlotA.YData = Taq;
            app.AnimationPlotA.XData = xq;

            title(app.UIAxes3, "t = " + num2str(app.AnimationIterations * period + t(app.CurrentFrame)))

            app.CurrentFrame = app.CurrentFrame + 1;

            drawnow limitrate
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: LEditField, NSpinner, OmegaEditField, 
        % ...and 5 other components
        function Update(app, event)
            N = app.NSpinner.Value;
            dt = app.TimeStepEditField.Value;
            L = app.LEditField.Value;
            omega = app.OmegaEditField.Value;
            tmax = app.tmaxEditField.Value;

            showAnalytical = app.ShowAnalyticalCheckBox.Value;
            showNumerical = app.ShowNumericalCheckBox.Value;
            showGrid = app.ShowGridCheckBox.Value;

            [x, t, T, Ta] = FE_Askisi_23(app, N, dt, L, omega, tmax);

            plot(app.UIAxes, 0, 0, "Visible", "off");
            plot(app.UIAxes2, 0, 0, "Visible", "off");

            disp(size(t))

            for i = 1 : 15 : length(t)

                if (showNumerical)
                    plot(app.UIAxes, x, T(:, i));
                end

                if (showAnalytical)
                    plot(app.UIAxes, x, Ta(:, i), "LineStyle", "--");
                end

                if (i == 1)
                    hold(app.UIAxes, "on");
                end
            end
            hold(app.UIAxes, "off");
            
            [X, Tgrid] = meshgrid(x, t);

            if (showNumerical)
                app.s1 = surf(app.UIAxes2, X, Tgrid, T');
                hold(app.UIAxes2, "on");

                app.s1.FaceAlpha = 0.7;
                app.s1.FaceColor = 'red';

                if (~showGrid)
                    app.s1.EdgeColor = "none";
                else
                    app.s1.EdgeColor = "black";
                end
            end

            if (showAnalytical)
                app.s2 = surf(app.UIAxes2, X, Tgrid, Ta');

                app.s2.FaceAlpha = 0.7;
                app.s2.FaceColor = 'blue';

                if (~showGrid)
                    app.s2.EdgeColor = "none";
                else
                    app.s2.EdgeColor = "black";
                end
            end
            hold(app.UIAxes2, "off");

            legend(app.UIAxes2, {"Numerical", "Analytical"})

            xlabel(app.UIAxes2, "x");
            ylabel(app.UIAxes2, "t");
            zlabel(app.UIAxes2, "T");

        end

        % Button pushed function: StartButton
        function StartAnimation(app, event)
            app.AnimationPlot = plot(app.UIAxes3, 0, 0);
            app.AnimationPlot.Visible = "off";

            hold(app.UIAxes3, "on");
            app.AnimationPlotA = plot(app.UIAxes3, 0, 0);
            app.AnimationPlotA.Visible = "off";
            hold(app.UIAxes3, "off");

            app.IsAnimating = true;

            app.AnimationTimer = timer( ...
                "ExecutionMode", "fixedRate", ...
                "Period", 0.05, ...
                "TimerFcn", @(~, ~) updateAnimation(app));

            ylim(app.UIAxes3, [-1, 1])
            % legend(app.UIAxes3, {"Numerical", "Analytical"})

            start(app.AnimationTimer);
        end

        % Button pushed function: StopButton
        function StopAnimation(app, event)
            if (isvalid(app.AnimationTimer))
                stop(app.AnimationTimer);
                delete(app.AnimationTimer);
            end
        end

        % Value changed function: InterpolationDropDown
        function UpdateInterpolation(app, event)
            app.Interpolation = app.InterpolationDropDown.Value;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1618 899];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.GridLayout);
            title(app.UIAxes2, 'Temperature Surface Plot')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.FontName = 'Garamond';
            app.UIAxes2.Box = 'on';
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.FontSize = 20;
            app.UIAxes2.Layout.Row = [10 17];
            app.UIAxes2.Layout.Column = [13 24];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Temperature Distribution')
            xlabel(app.UIAxes, 'x')
            ylabel(app.UIAxes, 'T')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Garamond';
            app.UIAxes.Box = 'on';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.FontSize = 20;
            app.UIAxes.Layout.Row = [2 9];
            app.UIAxes.Layout.Column = [10 21];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.GridLayout);
            title(app.UIAxes3, 'Temperature over time')
            xlabel(app.UIAxes3, 'x')
            ylabel(app.UIAxes3, 'T')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.FontName = 'Garamond';
            app.UIAxes3.Box = 'on';
            app.UIAxes3.XGrid = 'on';
            app.UIAxes3.YGrid = 'on';
            app.UIAxes3.FontSize = 15;
            app.UIAxes3.Layout.Row = [10 15];
            app.UIAxes3.Layout.Column = [2 8];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.GridLayout);
            app.UIAxes4.XTick = [];
            app.UIAxes4.YTick = [];
            app.UIAxes4.Box = 'on';
            app.UIAxes4.XGrid = 'on';
            app.UIAxes4.YGrid = 'on';
            app.UIAxes4.Layout.Row = [8 9];
            app.UIAxes4.Layout.Column = [2 8];

            % Create NSpinnerLabel
            app.NSpinnerLabel = uilabel(app.GridLayout);
            app.NSpinnerLabel.HorizontalAlignment = 'right';
            app.NSpinnerLabel.FontName = 'Garamond';
            app.NSpinnerLabel.FontSize = 20;
            app.NSpinnerLabel.FontWeight = 'bold';
            app.NSpinnerLabel.Layout.Row = 2;
            app.NSpinnerLabel.Layout.Column = 4;
            app.NSpinnerLabel.Text = 'N';

            % Create NSpinner
            app.NSpinner = uispinner(app.GridLayout);
            app.NSpinner.Limits = [2 Inf];
            app.NSpinner.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.NSpinner.HorizontalAlignment = 'center';
            app.NSpinner.FontName = 'Garamond';
            app.NSpinner.FontSize = 20;
            app.NSpinner.FontWeight = 'bold';
            app.NSpinner.Layout.Row = 2;
            app.NSpinner.Layout.Column = [5 6];
            app.NSpinner.Value = 2;

            % Create TimeStepEditFieldLabel
            app.TimeStepEditFieldLabel = uilabel(app.GridLayout);
            app.TimeStepEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeStepEditFieldLabel.FontName = 'Garamond';
            app.TimeStepEditFieldLabel.FontSize = 20;
            app.TimeStepEditFieldLabel.FontWeight = 'bold';
            app.TimeStepEditFieldLabel.Layout.Row = 3;
            app.TimeStepEditFieldLabel.Layout.Column = [3 4];
            app.TimeStepEditFieldLabel.Text = 'Time Step';

            % Create TimeStepEditField
            app.TimeStepEditField = uieditfield(app.GridLayout, 'numeric');
            app.TimeStepEditField.Limits = [0.0001 Inf];
            app.TimeStepEditField.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.TimeStepEditField.HorizontalAlignment = 'center';
            app.TimeStepEditField.FontName = 'Garamond';
            app.TimeStepEditField.FontSize = 20;
            app.TimeStepEditField.FontWeight = 'bold';
            app.TimeStepEditField.Layout.Row = 3;
            app.TimeStepEditField.Layout.Column = [5 6];
            app.TimeStepEditField.Value = 0.1;

            % Create LEditFieldLabel
            app.LEditFieldLabel = uilabel(app.GridLayout);
            app.LEditFieldLabel.HorizontalAlignment = 'right';
            app.LEditFieldLabel.FontName = 'Garamond';
            app.LEditFieldLabel.FontSize = 20;
            app.LEditFieldLabel.FontWeight = 'bold';
            app.LEditFieldLabel.Layout.Row = 4;
            app.LEditFieldLabel.Layout.Column = 4;
            app.LEditFieldLabel.Text = 'L';

            % Create LEditField
            app.LEditField = uieditfield(app.GridLayout, 'numeric');
            app.LEditField.Limits = [0.0001 Inf];
            app.LEditField.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.LEditField.HorizontalAlignment = 'center';
            app.LEditField.FontName = 'Garamond';
            app.LEditField.FontSize = 20;
            app.LEditField.FontWeight = 'bold';
            app.LEditField.Layout.Row = 4;
            app.LEditField.Layout.Column = [5 6];
            app.LEditField.Value = 5;

            % Create OmegaEditFieldLabel
            app.OmegaEditFieldLabel = uilabel(app.GridLayout);
            app.OmegaEditFieldLabel.HorizontalAlignment = 'right';
            app.OmegaEditFieldLabel.FontName = 'Garamond';
            app.OmegaEditFieldLabel.FontSize = 20;
            app.OmegaEditFieldLabel.FontWeight = 'bold';
            app.OmegaEditFieldLabel.Layout.Row = 5;
            app.OmegaEditFieldLabel.Layout.Column = [3 4];
            app.OmegaEditFieldLabel.Text = 'Omega';

            % Create OmegaEditField
            app.OmegaEditField = uieditfield(app.GridLayout, 'numeric');
            app.OmegaEditField.Limits = [0 Inf];
            app.OmegaEditField.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.OmegaEditField.HorizontalAlignment = 'center';
            app.OmegaEditField.FontName = 'Garamond';
            app.OmegaEditField.FontSize = 20;
            app.OmegaEditField.FontWeight = 'bold';
            app.OmegaEditField.Layout.Row = 5;
            app.OmegaEditField.Layout.Column = [5 6];
            app.OmegaEditField.Value = 2;

            % Create tmaxEditFieldLabel
            app.tmaxEditFieldLabel = uilabel(app.GridLayout);
            app.tmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.tmaxEditFieldLabel.FontName = 'Garamond';
            app.tmaxEditFieldLabel.FontSize = 20;
            app.tmaxEditFieldLabel.FontWeight = 'bold';
            app.tmaxEditFieldLabel.Layout.Row = 6;
            app.tmaxEditFieldLabel.Layout.Column = 4;
            app.tmaxEditFieldLabel.Text = 'tmax';

            % Create tmaxEditField
            app.tmaxEditField = uieditfield(app.GridLayout, 'numeric');
            app.tmaxEditField.Limits = [0 Inf];
            app.tmaxEditField.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.tmaxEditField.HorizontalAlignment = 'center';
            app.tmaxEditField.FontName = 'Garamond';
            app.tmaxEditField.FontSize = 20;
            app.tmaxEditField.FontWeight = 'bold';
            app.tmaxEditField.Layout.Row = 6;
            app.tmaxEditField.Layout.Column = [5 6];
            app.tmaxEditField.Value = 10;

            % Create ShowAnalyticalCheckBox
            app.ShowAnalyticalCheckBox = uicheckbox(app.GridLayout);
            app.ShowAnalyticalCheckBox.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.ShowAnalyticalCheckBox.Text = 'Show Analytical';
            app.ShowAnalyticalCheckBox.FontName = 'Garamond';
            app.ShowAnalyticalCheckBox.FontSize = 20;
            app.ShowAnalyticalCheckBox.FontWeight = 'bold';
            app.ShowAnalyticalCheckBox.Layout.Row = 7;
            app.ShowAnalyticalCheckBox.Layout.Column = [6 8];

            % Create ShowNumericalCheckBox
            app.ShowNumericalCheckBox = uicheckbox(app.GridLayout);
            app.ShowNumericalCheckBox.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.ShowNumericalCheckBox.Text = 'Show Numerical';
            app.ShowNumericalCheckBox.FontName = 'Garamond';
            app.ShowNumericalCheckBox.FontSize = 20;
            app.ShowNumericalCheckBox.FontWeight = 'bold';
            app.ShowNumericalCheckBox.Layout.Row = 7;
            app.ShowNumericalCheckBox.Layout.Column = [2 4];
            app.ShowNumericalCheckBox.Value = true;

            % Create ShowGridCheckBox
            app.ShowGridCheckBox = uicheckbox(app.GridLayout);
            app.ShowGridCheckBox.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.ShowGridCheckBox.Text = 'Show Grid';
            app.ShowGridCheckBox.FontName = 'Garamond';
            app.ShowGridCheckBox.FontSize = 20;
            app.ShowGridCheckBox.FontWeight = 'bold';
            app.ShowGridCheckBox.Layout.Row = 11;
            app.ShowGridCheckBox.Layout.Column = [25 27];
            app.ShowGridCheckBox.Value = true;

            % Create StartButton
            app.StartButton = uibutton(app.GridLayout, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartAnimation, true);
            app.StartButton.FontName = 'Garamond';
            app.StartButton.FontSize = 20;
            app.StartButton.FontWeight = 'bold';
            app.StartButton.Layout.Row = 16;
            app.StartButton.Layout.Column = [2 4];
            app.StartButton.Text = 'Start';

            % Create StopButton
            app.StopButton = uibutton(app.GridLayout, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopAnimation, true);
            app.StopButton.FontName = 'Garamond';
            app.StopButton.FontSize = 20;
            app.StopButton.FontWeight = 'bold';
            app.StopButton.Layout.Row = 16;
            app.StopButton.Layout.Column = [6 8];
            app.StopButton.Text = 'Stop';

            % Create InterpolationDropDownLabel
            app.InterpolationDropDownLabel = uilabel(app.GridLayout);
            app.InterpolationDropDownLabel.HorizontalAlignment = 'center';
            app.InterpolationDropDownLabel.FontName = 'Garamond';
            app.InterpolationDropDownLabel.FontSize = 15;
            app.InterpolationDropDownLabel.FontWeight = 'bold';
            app.InterpolationDropDownLabel.Layout.Row = 10;
            app.InterpolationDropDownLabel.Layout.Column = [9 10];
            app.InterpolationDropDownLabel.Text = 'Interpolation';

            % Create InterpolationDropDown
            app.InterpolationDropDown = uidropdown(app.GridLayout);
            app.InterpolationDropDown.Items = {'None', 'Linear', 'Nearest', 'Spline'};
            app.InterpolationDropDown.ItemsData = {'None', 'linear', 'nearest', 'spline'};
            app.InterpolationDropDown.ValueChangedFcn = createCallbackFcn(app, @UpdateInterpolation, true);
            app.InterpolationDropDown.FontName = 'Garamond';
            app.InterpolationDropDown.FontSize = 15;
            app.InterpolationDropDown.FontWeight = 'bold';
            app.InterpolationDropDown.Layout.Row = 10;
            app.InterpolationDropDown.Layout.Column = [11 12];
            app.InterpolationDropDown.Value = 'None';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Askisi_23_app_exported

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
