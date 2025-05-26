classdef Askisi_21_app_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        GridLayout                matlab.ui.container.GridLayout
        SecondValue               matlab.ui.control.NumericEditField
        ValueEditField_2Label     matlab.ui.control.Label
        SecondBCSelect            matlab.ui.control.DropDown
        ndBoundaryConditionx1DropDownLabel  matlab.ui.control.Label
        FirstValue                matlab.ui.control.NumericEditField
        ValueEditFieldLabel       matlab.ui.control.Label
        FirstBCSelect             matlab.ui.control.DropDown
        stBoundaryConditionLabel  matlab.ui.control.Label
        FD_ErrorLabel             matlab.ui.control.Label
        FE_ErrorLabel             matlab.ui.control.Label
        InitialLabel              matlab.ui.control.Label
        NSpinner                  matlab.ui.control.Spinner
        NSpinnerLabel             matlab.ui.control.Label
        LambdaEditField           matlab.ui.control.NumericEditField
        LambdaEditFieldLabel      matlab.ui.control.Label
        UIAxes                    matlab.ui.control.UIAxes
        Plot_2                    matlab.ui.control.UIAxes
        Plot_1                    matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        FirstBoundaryType % Description
        SecondBoundaryType
        FirstBoundaryValue
        SecondBoundaryValue
    end
    
    methods (Access = private)
        
        % =================== Finite Elements Method ====================
        function [Tn, x] = FEM_Askisi_21(app, N, lambda)
    
            Ne = N - 1;

            % Matrices
            A = zeros(N);
            Z = zeros(Ne, 2);

            % Vectors
            b = zeros(N, 1);
            beta = zeros(N, 1);
            x = zeros(N, 1);
            
            x(1) = 0;
            x(N) = 1;
            
            % Step size
            dx = 1. / double(N - 1);
            
            % Fill matrices
            for i = 2 : (N - 1)
                x(i) = (i - 1) * dx;

                % Beta vector
                beta(i) = exp(-lambda * x(i));
            end

            for ie = 1 : Ne
                Z(ie, 1) = ie;
                Z(ie, 2) = ie + 1;
            end

            Ll = 1 / dx * [
                1, -1 ;
                -1, 1
            ];

            Ml = dx / 6 * [
                2, 1 ;
                1, 2
            ];

            for ie = 1 : Ne
                for i = 1 : 2
                    il = Z(ie, i);

                    for j = 1 : 2
                        jl = Z(ie, j);

                        A(il, jl) = A(il, jl) + Ll(i, j);
                        b(il) = b(il) + Ml(i, j) * beta(jl);
                    end
                end
            end
            
            if (app.FirstBoundaryType == "Dirichlet")
                A(1, :) = 0.;
                A(1, 1) = 1.;
                b(1) = app.FirstBoundaryValue;
            else
                b(1) = b(1) - app.FirstBoundaryValue;
            end

            if (app.SecondBoundaryType == "Dirichlet")
                A(N, :) = 0.;
                A(N, N) = 1.;
                b(N) = app.SecondBoundaryValue;
            else
                b(N) = b(N) + app.SecondBoundaryValue;
            end

            Tn = A \ b;
        end





        % ================ Finite Differences Method ==================
        function [Tn, x] = FD_Askisi_21(app, N, lambda)
            
            % Matrix
            A = zeros(N);

            %Vectors
            x = zeros(N, 1);
            b = zeros(N, 1);

            dx = 1. / double(N - 1);

            for i = 1 : N
                x(i) = (i - 1) * dx;
            end

            % Fill matrices
            for i = 2 : (N - 1)
                A(i, i) = -2;
                A(i, i - 1) = 1;
                A(i, i + 1) = 1;

                b(i) = -dx^2 * exp(-lambda * x(i));
            end

            % Apply boundary conditions
            if (app.FirstBoundaryType == "Dirichlet")
                A(1, 1) = 1.;
            end

            if (app.FirstBoundaryType == "Neumann")
                A(1, 1) = -1.;
                A(1, 2) = 1.;
            end

            b(1) = app.FirstBoundaryValue;

            if (app.SecondBoundaryType == "Dirichlet")
                A(N, N) = 1.;
            end

            if (app.SecondBoundaryType == "Neumann")
                A(N, N - 1) = -1.;
                A(N, N) = 1.;
            end

            b(N) = app.SecondBoundaryValue;

            Tn = A \ b;
        end
        




        % ==================== Analytical Solution =====================
        function [Ta, x] = AN_Askisi_21(app, N, lambda)
            
            x = zeros(N, 1);

            dx = 1. / double(N - 1);

            for i = 1 : N
                x(i) = (i - 1) * dx;
            end

            % Solve for specified boundary conditions
            c1 = 0;
            c2 = 0;

            T1 = app.FirstBoundaryValue;
            T2 = app.SecondBoundaryValue;

            % Set integration constants accordingly
            if (app.FirstBoundaryType == "Neumann")
                if (app.SecondBoundaryType == "Dirichlet")
                    c1 = T1 - (1 / lambda);
                    c2 = T2 + (1 / lambda ^ 2) * exp(-lambda) - c1;
                end
            else
                if (app.SecondBoundaryType == "Dirichlet")
                    c1 = (T2 - T1) - (1 / lambda ^ 2) * (1 - exp(-lambda));
                    c2 = T1 + (1 / lambda ^ 2);
                else
                    c1 = T2 - (exp(-lambda) / lambda);
                    c2 = T1 + (1 / lambda ^ 2);
                end
            end

            % Calculate analytical solution
            Ta = -1 / (lambda ^ 2) * exp(-lambda * x) + c1 * x + c2;
        end





        % ======================== Update Plots =========================
        function [] = Update(app)

            % Get input data from application
            N = app.NSpinner.Value;
            lambda = app.LambdaEditField.Value;

            % Get boundary conditions
            % and make sure you can't select 2 Neumann
            if (app.FirstBCSelect.Value == "Neumann")
                app.SecondBCSelect.Value = "Dirichlet";
            end

            if (app.SecondBCSelect.Value == "Neumann")
                app.FirstBCSelect.Value = "Dirichlet";
            end

            app.FirstBoundaryType = app.FirstBCSelect.Value;
            app.FirstBoundaryValue = app.FirstValue.Value;

            app.SecondBoundaryType = app.SecondBCSelect.Value;
            app.SecondBoundaryValue = app.SecondValue.Value;

            % Calculate solutions
            [Ta, xa] = AN_Askisi_21(app, N, lambda);
            [Tn1, x1] = FEM_Askisi_21(app, N, lambda);
            [Tn2, x2] = FD_Askisi_21(app, N, lambda);

            app.InitialLabel.Visible = "off";

            % Display relative errors
            app.FE_ErrorLabel.Text = "Max relative error: " + num2str(max(abs(Ta - Tn1)) / max(Ta) * 100, 3) + "%";
            app.FD_ErrorLabel.Text = "Max relative error: " + num2str(max(abs(Ta - Tn2)) / max(Ta) * 100, 3) + "%";

            % Plot solutions
            plot(app.Plot_1, x1, Tn1);
            hold(app.Plot_1, "on")
            scatter(app.Plot_1, xa, Ta);
            hold(app.Plot_1, "off");

            MaxVal1 = max(max(Ta), max(Tn1));
            app.Plot_1.YLim = [0, MaxVal1];
            app.Plot_1.YTick = 0:(MaxVal1 / 10):MaxVal1;

            plot(app.Plot_2, x2, Tn2);
            hold(app.Plot_2, "on")
            scatter(app.Plot_2, xa, Ta);
            hold(app.Plot_2, "off");

            MaxVal2 = max(max(Ta), max(Tn2));
            app.Plot_2.YLim = [0, MaxVal2];
            app.Plot_2.YTick = 0:(MaxVal2 / 10):MaxVal2;

            title(app.Plot_1, "Finite Elements (N = " + num2str(N) + ")");
            title(app.Plot_2, "Finite Differences (N = " + num2str(N) + ")");

            legend(app.Plot_1, {"Numerical", "Analytical"}, "Location", "southeast")
            legend(app.Plot_2, {"Numerical", "Analytical"}, "Location", "southeast")

            % ================== Update error graph ===================
            Nvals = 10:10:100;

            % Error vectors
            fd_errors = zeros(length(Nvals), 1);
            fe_errors = zeros(length(Nvals), 1);

            % Set up axis limits
            app.UIAxes.XLim = [0, 100];

            % FIll error vectors
            for i = 1 : length(Nvals)
                Ta = app.AN_Askisi_21(Nvals(i), lambda);
                T1 = app.FD_Askisi_21(Nvals(i), lambda);
                T2 = app.FEM_Askisi_21(Nvals(i), lambda);

                fd_errors(i) = max(abs(T1 - Ta) / max(Ta));
                fe_errors(i) = max(abs(T2 - Ta) / max(Ta));
            end

            % Plot graph (logarithmic)
            loglog(app.UIAxes, Nvals, fd_errors);
            hold(app.UIAxes, "on");
            loglog(app.UIAxes, Nvals, fe_errors);
            hold(app.UIAxes, "off");

            legend(app.UIAxes, {"Finite Differences", "Finite Elements"})
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: FirstBCSelect, FirstValue, 
        % ...and 4 other components
        function LambdaValueChanged(app, event)
            app.Update();
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1623 930];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};

            % Create Plot_1
            app.Plot_1 = uiaxes(app.GridLayout);
            title(app.Plot_1, 'Finite Elements')
            xlabel(app.Plot_1, 'x')
            ylabel(app.Plot_1, 'Temperature')
            app.Plot_1.FontName = 'Garamond';
            app.Plot_1.XTick = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
            app.Plot_1.YTick = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
            app.Plot_1.Box = 'on';
            app.Plot_1.XGrid = 'on';
            app.Plot_1.YGrid = 'on';
            app.Plot_1.FontSize = 20;
            app.Plot_1.Layout.Row = [2 7];
            app.Plot_1.Layout.Column = [2 7];

            % Create Plot_2
            app.Plot_2 = uiaxes(app.GridLayout);
            title(app.Plot_2, 'Finite Differences')
            xlabel(app.Plot_2, 'x')
            ylabel(app.Plot_2, 'Temperature')
            app.Plot_2.FontName = 'Garamond';
            app.Plot_2.XTick = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
            app.Plot_2.YTick = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
            app.Plot_2.Box = 'on';
            app.Plot_2.XGrid = 'on';
            app.Plot_2.YGrid = 'on';
            app.Plot_2.FontSize = 20;
            app.Plot_2.Layout.Row = [2 7];
            app.Plot_2.Layout.Column = [10 15];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Error - N')
            xlabel(app.UIAxes, 'N')
            ylabel(app.UIAxes, 'Max relative error')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Garamond';
            app.UIAxes.Box = 'on';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.FontSize = 20;
            app.UIAxes.Layout.Row = [9 14];
            app.UIAxes.Layout.Column = [6 12];

            % Create LambdaEditFieldLabel
            app.LambdaEditFieldLabel = uilabel(app.GridLayout);
            app.LambdaEditFieldLabel.HorizontalAlignment = 'center';
            app.LambdaEditFieldLabel.FontName = 'Garamond';
            app.LambdaEditFieldLabel.FontSize = 20;
            app.LambdaEditFieldLabel.FontWeight = 'bold';
            app.LambdaEditFieldLabel.Layout.Row = 5;
            app.LambdaEditFieldLabel.Layout.Column = [8 9];
            app.LambdaEditFieldLabel.Text = 'Lambda';

            % Create LambdaEditField
            app.LambdaEditField = uieditfield(app.GridLayout, 'numeric');
            app.LambdaEditField.Limits = [0.001 Inf];
            app.LambdaEditField.ValueChangedFcn = createCallbackFcn(app, @LambdaValueChanged, true);
            app.LambdaEditField.HorizontalAlignment = 'center';
            app.LambdaEditField.FontName = 'Garamond';
            app.LambdaEditField.FontSize = 20;
            app.LambdaEditField.FontWeight = 'bold';
            app.LambdaEditField.Layout.Row = 6;
            app.LambdaEditField.Layout.Column = [8 9];
            app.LambdaEditField.Value = 1;

            % Create NSpinnerLabel
            app.NSpinnerLabel = uilabel(app.GridLayout);
            app.NSpinnerLabel.HorizontalAlignment = 'center';
            app.NSpinnerLabel.FontName = 'Garamond';
            app.NSpinnerLabel.FontSize = 20;
            app.NSpinnerLabel.FontWeight = 'bold';
            app.NSpinnerLabel.Layout.Row = 2;
            app.NSpinnerLabel.Layout.Column = [8 9];
            app.NSpinnerLabel.Text = 'N';

            % Create NSpinner
            app.NSpinner = uispinner(app.GridLayout);
            app.NSpinner.Limits = [2 Inf];
            app.NSpinner.ValueChangedFcn = createCallbackFcn(app, @LambdaValueChanged, true);
            app.NSpinner.HorizontalAlignment = 'center';
            app.NSpinner.FontName = 'Garamond';
            app.NSpinner.FontSize = 20;
            app.NSpinner.FontWeight = 'bold';
            app.NSpinner.Layout.Row = 3;
            app.NSpinner.Layout.Column = [8 9];
            app.NSpinner.Value = 2;

            % Create InitialLabel
            app.InitialLabel = uilabel(app.GridLayout);
            app.InitialLabel.HorizontalAlignment = 'center';
            app.InitialLabel.FontName = 'Garamond';
            app.InitialLabel.FontSize = 20;
            app.InitialLabel.Layout.Row = 8;
            app.InitialLabel.Layout.Column = [5 12];
            app.InitialLabel.Text = 'Modify any value to plot.';

            % Create FE_ErrorLabel
            app.FE_ErrorLabel = uilabel(app.GridLayout);
            app.FE_ErrorLabel.HorizontalAlignment = 'center';
            app.FE_ErrorLabel.VerticalAlignment = 'top';
            app.FE_ErrorLabel.FontName = 'Garamond';
            app.FE_ErrorLabel.FontSize = 20;
            app.FE_ErrorLabel.Layout.Row = 8;
            app.FE_ErrorLabel.Layout.Column = [3 6];
            app.FE_ErrorLabel.Text = '';

            % Create FD_ErrorLabel
            app.FD_ErrorLabel = uilabel(app.GridLayout);
            app.FD_ErrorLabel.HorizontalAlignment = 'center';
            app.FD_ErrorLabel.VerticalAlignment = 'top';
            app.FD_ErrorLabel.FontName = 'Garamond';
            app.FD_ErrorLabel.FontSize = 20;
            app.FD_ErrorLabel.Layout.Row = 8;
            app.FD_ErrorLabel.Layout.Column = [11 14];
            app.FD_ErrorLabel.Text = '';

            % Create stBoundaryConditionLabel
            app.stBoundaryConditionLabel = uilabel(app.GridLayout);
            app.stBoundaryConditionLabel.HorizontalAlignment = 'center';
            app.stBoundaryConditionLabel.FontName = 'Garamond';
            app.stBoundaryConditionLabel.FontSize = 20;
            app.stBoundaryConditionLabel.FontWeight = 'bold';
            app.stBoundaryConditionLabel.Layout.Row = 9;
            app.stBoundaryConditionLabel.Layout.Column = [1 2];
            app.stBoundaryConditionLabel.Text = {'1st Boundary'; 'Condition (x = 0)'};

            % Create FirstBCSelect
            app.FirstBCSelect = uidropdown(app.GridLayout);
            app.FirstBCSelect.Items = {'Dirichlet', 'Neumann'};
            app.FirstBCSelect.ValueChangedFcn = createCallbackFcn(app, @LambdaValueChanged, true);
            app.FirstBCSelect.FontName = 'Garamond';
            app.FirstBCSelect.FontSize = 20;
            app.FirstBCSelect.Layout.Row = 9;
            app.FirstBCSelect.Layout.Column = [3 4];
            app.FirstBCSelect.Value = 'Dirichlet';

            % Create ValueEditFieldLabel
            app.ValueEditFieldLabel = uilabel(app.GridLayout);
            app.ValueEditFieldLabel.HorizontalAlignment = 'center';
            app.ValueEditFieldLabel.FontName = 'Garamond';
            app.ValueEditFieldLabel.FontSize = 20;
            app.ValueEditFieldLabel.FontWeight = 'bold';
            app.ValueEditFieldLabel.Layout.Row = 10;
            app.ValueEditFieldLabel.Layout.Column = [1 2];
            app.ValueEditFieldLabel.Text = 'Value';

            % Create FirstValue
            app.FirstValue = uieditfield(app.GridLayout, 'numeric');
            app.FirstValue.ValueChangedFcn = createCallbackFcn(app, @LambdaValueChanged, true);
            app.FirstValue.HorizontalAlignment = 'center';
            app.FirstValue.FontName = 'Garamond';
            app.FirstValue.FontSize = 20;
            app.FirstValue.FontWeight = 'bold';
            app.FirstValue.Layout.Row = 10;
            app.FirstValue.Layout.Column = [3 4];

            % Create ndBoundaryConditionx1DropDownLabel
            app.ndBoundaryConditionx1DropDownLabel = uilabel(app.GridLayout);
            app.ndBoundaryConditionx1DropDownLabel.HorizontalAlignment = 'center';
            app.ndBoundaryConditionx1DropDownLabel.FontName = 'Garamond';
            app.ndBoundaryConditionx1DropDownLabel.FontSize = 20;
            app.ndBoundaryConditionx1DropDownLabel.FontWeight = 'bold';
            app.ndBoundaryConditionx1DropDownLabel.Layout.Row = 12;
            app.ndBoundaryConditionx1DropDownLabel.Layout.Column = [1 2];
            app.ndBoundaryConditionx1DropDownLabel.Text = {'2nd Boundary'; 'Condition (x = 1)'};

            % Create SecondBCSelect
            app.SecondBCSelect = uidropdown(app.GridLayout);
            app.SecondBCSelect.Items = {'Dirichlet', 'Neumann'};
            app.SecondBCSelect.ValueChangedFcn = createCallbackFcn(app, @LambdaValueChanged, true);
            app.SecondBCSelect.FontName = 'Garamond';
            app.SecondBCSelect.FontSize = 20;
            app.SecondBCSelect.Layout.Row = 12;
            app.SecondBCSelect.Layout.Column = [3 4];
            app.SecondBCSelect.Value = 'Dirichlet';

            % Create ValueEditField_2Label
            app.ValueEditField_2Label = uilabel(app.GridLayout);
            app.ValueEditField_2Label.HorizontalAlignment = 'center';
            app.ValueEditField_2Label.FontName = 'Garamond';
            app.ValueEditField_2Label.FontSize = 20;
            app.ValueEditField_2Label.FontWeight = 'bold';
            app.ValueEditField_2Label.Layout.Row = 13;
            app.ValueEditField_2Label.Layout.Column = [1 2];
            app.ValueEditField_2Label.Text = 'Value';

            % Create SecondValue
            app.SecondValue = uieditfield(app.GridLayout, 'numeric');
            app.SecondValue.ValueChangedFcn = createCallbackFcn(app, @LambdaValueChanged, true);
            app.SecondValue.HorizontalAlignment = 'center';
            app.SecondValue.FontName = 'Garamond';
            app.SecondValue.FontSize = 20;
            app.SecondValue.FontWeight = 'bold';
            app.SecondValue.Layout.Row = 13;
            app.SecondValue.Layout.Column = [3 4];
            app.SecondValue.Value = 1;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Askisi_21_app_exported

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
