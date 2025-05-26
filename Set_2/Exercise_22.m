classdef Askisi_22_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        GridLayout              matlab.ui.container.GridLayout
        TimeStepEditField       matlab.ui.control.NumericEditField
        TimeStepEditFieldLabel  matlab.ui.control.Label
        LambdaEditField         matlab.ui.control.NumericEditField
        LambdaEditFieldLabel    matlab.ui.control.Label
        UIAxes3                 matlab.ui.control.UIAxes
        UIAxes2                 matlab.ui.control.UIAxes
        UIAxes                  matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        T % Description
        tmax
    end
    
    methods (Access = public)
        
        function [x, y] = ImpicitEuler(~, lambda, dt, tmax)
            
            N = round(tmax / dt);

            A = zeros(2, 2);
            f = zeros(2, N);
            
            A(1, 1) = 1;
            A(1, 2) = dt;
            A(2, 1) = -lambda^2 * dt;
            A(2, 2) = 1;
            
            f(1, 1) = 0;
            f(2, 1) = lambda;
            
            for i = 2 : N
                f(:, i) = A * f(:, i - 1);
            end

            y = f(1, :);
            x = linspace(0, tmax, N);
        end

        function [x, y] = Analytical(~, lambda, dt, tmax)

            N = round(tmax / dt);

            x = linspace(0, tmax, N);
            y = sin(lambda * x);
        end
    end
    
    methods (Access = private)
        
        function [x, y] = ExplicitEuler(~, lambda, dt, tmax)
            
            N = round(tmax / dt);

            A = zeros(2, 2);
            f = zeros(2, N);
            
            A(1, 1) = 1;
            A(1, 2) = dt;
            A(2, 1) = -lambda^2 * dt;
            A(2, 2) = 1;
            
            A = A / (1 + (lambda^2 * dt^2));

            f(1, 1) = 0;
            f(2, 1) = lambda;
            
            for i = 2 : N
                f(:, i) = A * f(:, i - 1);
            end

            y = f(1, :);
            x = linspace(0, tmax, N);
        end
        
        function [x, y] = CrankNicholson(~, lambda, dt, tmax)

            N = round(tmax / dt);

            A = zeros(2, 2);
            f = zeros(2, N);
            
            A(1, 1) = 2;
            A(1, 2) = dt;
            A(2, 1) = -lambda^2 * dt;
            A(2, 2) = 2;
            
            A = A^2;

            A = A / (4 + (lambda^2 * dt^2));

            f(1, 1) = 0;
            f(2, 1) = lambda;
            
            for i = 2 : N
                f(:, i) = A * f(:, i - 1);
            end

            y = f(1, :);
            x = linspace(0, tmax, N);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: LambdaEditField, TimeStepEditField
        function Update(app, event)
            lambda = app.LambdaEditField.Value;
            dt = app.TimeStepEditField.Value;

            app.T = 2 * pi / lambda;
            app.tmax = 5 * app.T;

            [x_ie, y_ie] = ImpicitEuler(app, lambda, dt, app.tmax);
            [x_ee, y_ee] = ExplicitEuler(app, lambda, dt, app.tmax);
            [x_cn, y_cn] = CrankNicholson(app, lambda, dt, app.tmax);
            [x_a, y_a] = Analytical(app, lambda, dt, app.tmax);

            app.UIAxes.YLim = [min(y_ie), max(y_ie)];
            plot(app.UIAxes, x_ie, y_ie);
            hold(app.UIAxes, "on");
            plot(app.UIAxes, x_a, y_a, "LineStyle", "--");
            hold(app.UIAxes, "off");

            title(app.UIAxes, "Implicit Euler (Lambda = " + num2str(lambda) + ", Time Step = " + num2str(dt) + ")");
            legend(app.UIAxes, {"Implicit Euler", "Analytical"}, "Location", "southwest");

            app.UIAxes2.YLim = [min(min(y_ee), min(y_a)), max(max(y_ee), max(y_a))];
            plot(app.UIAxes2, x_ee, y_ee);
            hold(app.UIAxes2, "on");
            plot(app.UIAxes2, x_a, y_a, "LineStyle", "--");
            hold(app.UIAxes2, "off");

            title(app.UIAxes2, "Explicit Euler (Lambda = " + num2str(lambda) + ", Time Step = " + num2str(dt) + ")");
            legend(app.UIAxes2, {"Explicit Euler", "Analytical"}, "Location", "southwest");

            app.UIAxes3.YLim = [-1, 1];
            plot(app.UIAxes3, x_cn, y_cn);
            hold(app.UIAxes3, "on");
            plot(app.UIAxes3, x_a, y_a, "LineStyle", "--");
            hold(app.UIAxes3, "off");

            title(app.UIAxes3, "Crank - Nicholson (Lambda = " + num2str(lambda) + ", Time Step = " + num2str(dt) + ")");
            legend(app.UIAxes3, {"Crank - Nicholson", "Analytical"}, "Location", "southwest");
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1577 859];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Implicit Euler')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Garamond';
            app.UIAxes.Box = 'on';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.FontSize = 20;
            app.UIAxes.Layout.Row = [2 6];
            app.UIAxes.Layout.Column = [2 7];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.GridLayout);
            title(app.UIAxes2, 'Explicit Euler')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.FontName = 'garamond';
            app.UIAxes2.Box = 'on';
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.FontSize = 20;
            app.UIAxes2.Layout.Row = [2 6];
            app.UIAxes2.Layout.Column = [12 17];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.GridLayout);
            title(app.UIAxes3, 'Crank - Nicholson')
            xlabel(app.UIAxes3, 'X')
            ylabel(app.UIAxes3, 'Y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.FontName = 'Garamond';
            app.UIAxes3.Box = 'on';
            app.UIAxes3.XGrid = 'on';
            app.UIAxes3.YGrid = 'on';
            app.UIAxes3.FontSize = 20;
            app.UIAxes3.Layout.Row = [7 11];
            app.UIAxes3.Layout.Column = [7 12];

            % Create LambdaEditFieldLabel
            app.LambdaEditFieldLabel = uilabel(app.GridLayout);
            app.LambdaEditFieldLabel.HorizontalAlignment = 'center';
            app.LambdaEditFieldLabel.FontName = 'Garamond';
            app.LambdaEditFieldLabel.FontSize = 20;
            app.LambdaEditFieldLabel.FontWeight = 'bold';
            app.LambdaEditFieldLabel.Layout.Row = 3;
            app.LambdaEditFieldLabel.Layout.Column = 9;
            app.LambdaEditFieldLabel.Text = 'Lambda';

            % Create LambdaEditField
            app.LambdaEditField = uieditfield(app.GridLayout, 'numeric');
            app.LambdaEditField.Limits = [0.0001 Inf];
            app.LambdaEditField.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.LambdaEditField.HorizontalAlignment = 'center';
            app.LambdaEditField.FontName = 'Garamond';
            app.LambdaEditField.FontSize = 20;
            app.LambdaEditField.FontWeight = 'bold';
            app.LambdaEditField.Layout.Row = 3;
            app.LambdaEditField.Layout.Column = 10;
            app.LambdaEditField.Value = 2;

            % Create TimeStepEditFieldLabel
            app.TimeStepEditFieldLabel = uilabel(app.GridLayout);
            app.TimeStepEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeStepEditFieldLabel.FontName = 'Garamond';
            app.TimeStepEditFieldLabel.FontSize = 20;
            app.TimeStepEditFieldLabel.FontWeight = 'bold';
            app.TimeStepEditFieldLabel.Layout.Row = 4;
            app.TimeStepEditFieldLabel.Layout.Column = [8 9];
            app.TimeStepEditFieldLabel.Text = 'Time Step';

            % Create TimeStepEditField
            app.TimeStepEditField = uieditfield(app.GridLayout, 'numeric');
            app.TimeStepEditField.Limits = [0.0001 1];
            app.TimeStepEditField.ValueChangedFcn = createCallbackFcn(app, @Update, true);
            app.TimeStepEditField.HorizontalAlignment = 'center';
            app.TimeStepEditField.FontName = 'Garamond';
            app.TimeStepEditField.FontSize = 20;
            app.TimeStepEditField.FontWeight = 'bold';
            app.TimeStepEditField.Layout.Row = 4;
            app.TimeStepEditField.Layout.Column = 10;
            app.TimeStepEditField.Value = 1;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Askisi_22_app_exported

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
