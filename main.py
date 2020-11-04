import numpy as np
import math
from PyQt5 import QtWidgets, uic
import pyqtgraph as pg
import sys

app = QtWidgets.QApplication(sys.argv)


class DiffEquation:

    def __init__(self, window):
        self.x_0 = window.x_0.value()
        self.y_0 = window.y_0.value()
        self.X = window.X.value()
        self.N = int(window.N.value())
        self.h = (self.X - self.x_0) / self.N
        self.N_start = int(window.N_start.value())
        self.N_finish = int(window.N_finish.value())

    def update_values(self, x_0=None, y_0=None, X=None, N=None):
        if x_0:
            self.x_0 = x_0
            self.h = (self.X - self.x_0) / self.N
        if y_0:
            self.y_0 = y_0
        if X:
            self.X = X
            self.h = (self.X - self.x_0) / self.N
        if N:
            self.N = N
            self.h = (self.X - self.x_0) / self.N

    def get_list_of_x(self):
        return np.arange(self.x_0, self.X + self.h/2, self.h)

    def get_n_list_of_error(self):
        return np.arange(self.N_start, self.N_finish + 1, 1)

    def get_list_of_y(self, list_of_x):
        pass

    def get_list_of_global_errors(self, window):
        pass

    def get_list_of_local_errors(self, window):
        pass

    def get_list_of_max_global_errors(self, window):
        pass

    def get_list_of_max_local_errors(self, window):
        pass

    def find_f(self, x, y):
        pass


class ExactMethod(DiffEquation):

    def get_C(self):
        return (self.y_0 + math.log(self.x_0, math.e)) / math.pow(math.log(self.x_0, math.e), 2)

    def get_list_of_y(self, list_of_x):
        return [math.log(x, math.e) * (self.get_C() * math.log(x, math.e) - 1) for x in list_of_x]


class EulerMethod(DiffEquation):

    def find_f(self, x, y):
        return 1 / x + (2 * y) / (x * math.log(x, math.e))

    def get_list_of_y(self, list_of_x):
        y = self.y_0
        list_of_y = [y]
        for i in range(1, len(list_of_x)):
            y_i = list_of_y[i - 1] + self.h * self.find_f(list_of_x[i - 1], list_of_y[i - 1])
            list_of_y.append(y_i)
        return list_of_y

    def get_list_of_global_errors(self, window):
        exact = ExactMethod(window)
        euler = EulerMethod(window)
        exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
        euler_list_of_y = euler.get_list_of_y(exact.get_list_of_x())
        error_list = []
        for i in range(0, len(exact.get_list_of_x())):
            error_list.append(abs(exact_list_of_y[i] - euler_list_of_y[i]))
        return error_list

    def get_list_of_local_errors(self, window):
        exact = ExactMethod(window)
        exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
        error_list = [0.0]
        for i in range(0, len(exact.get_list_of_x()) - 1):
            x_i = self.x_0 + self.h * i
            error_list.append(abs(exact_list_of_y[i+1] - exact_list_of_y[i] - self.h *
                                  self.find_f(x_i, exact_list_of_y[i])))
        return error_list

    def get_list_of_max_global_errors(self, window):
        list_of_max_global_errors = []
        exact = ExactMethod(window)
        euler = EulerMethod(window)
        for j in self.get_n_list_of_error():
            exact.update_values(N=j)
            euler.update_values(N=j)
            exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
            euler_list_of_y = euler.get_list_of_y(exact.get_list_of_x())
            error_list = []
            for i in range(0, len(exact.get_list_of_x())):
                error_list.append(abs(exact_list_of_y[i] - euler_list_of_y[i]))
            list_of_max_global_errors.append(max(error_list))
        return list_of_max_global_errors

    def get_list_of_max_local_errors(self, window):
        list_of_max_local_errors = []
        exact = ExactMethod(window)
        for j in self.get_n_list_of_error():
            exact.update_values(N=j)
            exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
            error_list = [0.0]
            for i in range(0, len(exact.get_list_of_x()) - 1):
                x_i = exact.x_0 + exact.h * i
                error_list.append(
                    abs(exact_list_of_y[i + 1] - exact_list_of_y[i] - exact.h * self.find_f(x_i, exact_list_of_y[i])))
            list_of_max_local_errors.append(max(error_list))
        return list_of_max_local_errors


class ImprovedEulerMethod(DiffEquation):

    def find_f(self, x, y):
        return 1 / x + (2 * y) / (x * math.log(x, math.e))

    def get_list_of_y(self, list_of_x):
        y = self.y_0
        list_of_y = [y]
        for i in range(1, len(list_of_x)):
            K1 = self.find_f(list_of_x[i - 1], list_of_y[i - 1])
            K2 = self.find_f(list_of_x[i - 1] + self.h, list_of_y[i - 1] + self.h * K1)
            y_i = list_of_y[i - 1] + self.h * (K1 + K2) / 2
            list_of_y.append(y_i)
        return list_of_y

    def get_list_of_global_errors(self, window):
        exact = ExactMethod(window)
        improved_euler = ImprovedEulerMethod(window)
        exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
        improved_euler_list_of_y = improved_euler.get_list_of_y(exact.get_list_of_x())
        error_list = []
        for i in range(0, len(exact.get_list_of_x())):
            error_list.append(abs(exact_list_of_y[i] - improved_euler_list_of_y[i]))
        return error_list

    def get_list_of_local_errors(self, window):
        exact = ExactMethod(window)
        exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
        error_list = [0.0]
        for i in range(0, len(exact.get_list_of_x()) - 1):
            x_i = self.x_0 + self.h * i
            error_list.append(abs(exact_list_of_y[i + 1] - (exact_list_of_y[i] + self.h * (
                        self.find_f(x_i, exact_list_of_y[i]) + self.find_f(x_i + self.h,
                                                                           exact_list_of_y[i] + self.h * self.find_f(
                                                                               x_i, exact_list_of_y[i])))/2)))

        return error_list

    def get_list_of_max_global_errors(self, window):
        list_of_max_global_errors = []
        exact = ExactMethod(window)
        improved_euler = ImprovedEulerMethod(window)
        for j in self.get_n_list_of_error():
            exact.update_values(N=j)
            improved_euler.update_values(N=j)
            exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
            improved_euler_list_of_y = improved_euler.get_list_of_y(exact.get_list_of_x())
            error_list = []
            for i in range(0, len(exact.get_list_of_x())):
                error_list.append(abs(exact_list_of_y[i] - improved_euler_list_of_y[i]))
            list_of_max_global_errors.append(max(error_list))
        return list_of_max_global_errors

    def get_list_of_max_local_errors(self, window):
        list_of_max_local_errors = []
        exact = ExactMethod(window)
        for j in self.get_n_list_of_error():
            exact.update_values(N=j)
            exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
            error_list = [0.0]
            for i in range(0, len(exact.get_list_of_x()) - 1):
                x_i = exact.x_0 + exact.h * i
                error_list.append(abs(exact_list_of_y[i + 1] - (exact_list_of_y[i] + exact.h * (
                        self.find_f(x_i, exact_list_of_y[i]) + self.find_f(x_i + exact.h,
                                                                           exact_list_of_y[i] + exact.h * self.find_f(
                                                                               x_i, exact_list_of_y[i]))) / 2)))
            list_of_max_local_errors.append(max(error_list))
        return list_of_max_local_errors


class RungeKuttaMethod(DiffEquation):

    def find_f(self, x, y):
        return 1 / x + (2 * y) / (x * math.log(x, math.e))

    def get_list_of_y(self, list_of_x):
        y = self.y_0
        list_of_y = [y]
        for i in range(1, len(list_of_x)):
            K1 = self.find_f(list_of_x[i - 1], list_of_y[i - 1])
            K2 = self.find_f(list_of_x[i - 1] + self.h / 2, list_of_y[i - 1] + self.h * K1 / 2)
            K3 = self.find_f(list_of_x[i - 1] + self.h / 2, list_of_y[i - 1] + self.h * K2 / 2)
            K4 = self.find_f(list_of_x[i - 1] + self.h, list_of_y[i - 1] + self.h * K3)
            y_i = list_of_y[i - 1] + self.h * (K1 + 2 * K2 + 2 * K3 + K4) / 6
            list_of_y.append(y_i)
        return list_of_y

    def get_list_of_global_errors(self, window):
        exact = ExactMethod(window)
        runge_kutta = RungeKuttaMethod(window)
        exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
        runge_kutta_list_of_y = runge_kutta.get_list_of_y(exact.get_list_of_x())
        error_list = []
        for i in range(0, len(exact.get_list_of_x())):
            error_list.append(abs(exact_list_of_y[i] - runge_kutta_list_of_y[i]))
        return error_list

    def get_list_of_local_errors(self, window):
        exact = ExactMethod(window)
        exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
        error_list = [0.0]
        for i in range(0, len(exact.get_list_of_x()) - 1):
            x_i = self.x_0 + self.h * i
            K1 = self.find_f(x_i, exact_list_of_y[i])
            K2 = self.find_f(x_i + self.h / 2, exact_list_of_y[i] + self.h * K1 / 2)
            K3 = self.find_f(x_i + self.h / 2, exact_list_of_y[i] + self.h * K2 / 2)
            K4 = self.find_f(x_i + self.h, exact_list_of_y[i] + self.h * K3)
            y_i_plus_1 = exact_list_of_y[i] + self.h * (K1 + 2 * K2 + 2 * K3 + K4) / 6
            error_list.append(abs(exact_list_of_y[i + 1] - y_i_plus_1))
        return error_list

    def get_list_of_max_global_errors(self, window):
        list_of_max_global_errors = []
        exact = ExactMethod(window)
        runge_kutta = RungeKuttaMethod(window)
        for j in self.get_n_list_of_error():
            exact.update_values(N=j)
            runge_kutta.update_values(N=j)
            exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
            runge_kutta_list_of_y = runge_kutta.get_list_of_y(exact.get_list_of_x())
            error_list = []
            for i in range(0, len(exact.get_list_of_x())):
                error_list.append(abs(exact_list_of_y[i] - runge_kutta_list_of_y[i]))
            list_of_max_global_errors.append(max(error_list))
        return list_of_max_global_errors

    def get_list_of_max_local_errors(self, window):
        list_of_max_local_errors = []
        exact = ExactMethod(window)
        for j in self.get_n_list_of_error():
            exact.update_values(N=j)
            exact_list_of_y = exact.get_list_of_y(exact.get_list_of_x())
            error_list = [0.0]
            for i in range(0, len(exact.get_list_of_x()) - 1):
                x_i = exact.x_0 + exact.h * i
                K1 = self.find_f(x_i, exact_list_of_y[i])
                K2 = self.find_f(x_i + exact.h / 2, exact_list_of_y[i] + exact.h * K1 / 2)
                K3 = self.find_f(x_i + exact.h / 2, exact_list_of_y[i] + exact.h * K2 / 2)
                K4 = self.find_f(x_i + exact.h, exact_list_of_y[i] + exact.h * K3)
                y_i_plus_1 = exact_list_of_y[i] + exact.h * (K1 + 2 * K2 + 2 * K3 + K4) / 6
                error_list.append(abs(exact_list_of_y[i + 1] - y_i_plus_1))
            list_of_max_local_errors.append(max(error_list))
        return list_of_max_local_errors


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        uic.loadUi('./form.ui', self)

        # connecting all buttons, checkboxes and radio buttons
        self.pushButton.clicked.connect(self.plot_new_graphs)

        self.exact_checkBox.stateChanged.connect(self.show_or_hide_curve)
        self.euler_checkBox.stateChanged.connect(self.show_or_hide_curve)
        self.improved_euler_checkBox.stateChanged.connect(self.show_or_hide_curve)
        self.runge_kutta_checkBox.stateChanged.connect(self.show_or_hide_curve)

        self.LTE_button.toggled.connect(self.show_or_hide_errors_curve)
        self.GTE_button.toggled.connect(self.show_or_hide_errors_curve)

        # assigning titles and labels to graphs
        self.Solutions.setTitle("Solutions", color='k', size='30pt')
        self.Solutions.setLabel('left', 'y', color='k', **{'font-size': '25pt'})
        self.Solutions.setLabel('bottom', 'x', color='k', **{'font-size': '25pt'})

        self.Errors.setTitle("Errors", color='000', size='25pt')
        self.Errors.setLabel('left', 'Error', color='000', **{'font-size': '25pt'})
        self.Errors.setLabel('bottom', 'x', color='000', **{'font-size': '25pt'})

        self.Max_Errors.setTitle("Max Errors", color='000', size='25pt')
        self.Max_Errors.setLabel('left', 'Max Error', color='000', **{'font-size': '25pt'})
        self.Max_Errors.setLabel('bottom', 'N', color='000', **{'font-size': '25pt'})

        self.plot_new_graphs()

    # to set the parameters of the graphs and curves and plot them
    def plot(self, graph, hour=None, temperature=None, color=None, name=None):
        if hour is None or temperature is None:
            hour, temperature, color = [], [], (1, 0, 0)

        graph.setBackground('w')
        graph.showGrid(x=True, y=True)

        graph.plot(hour, temperature, pen=pg.mkPen(color=color, width=5), name=name)

    def plot_new_graphs(self):

        values = DiffEquation(self)

        # clear all graphs
        self.clear_graphs()

        self.Solutions.addLegend()
        self.Solutions.enableAutoRange()
        self.Errors.addLegend()
        self.Errors.enableAutoRange()
        self.Max_Errors.addLegend()
        self.Max_Errors.enableAutoRange()

        # plot solution for ExactMethod
        exact_solution = ExactMethod(self)
        self.plot(self.Solutions, exact_solution.get_list_of_x(),
                  exact_solution.get_list_of_y(exact_solution.get_list_of_x()), color=(255, 20, 147), name="exact")

        # plot solution, GTE and Max GTE curves for EulerMethod
        euler_solution = EulerMethod(self)
        self.plot(self.Solutions, euler_solution.get_list_of_x(),
                  euler_solution.get_list_of_y(euler_solution.get_list_of_x()), color=(255, 165, 0), name="euler")
        self.plot(self.Errors, euler_solution.get_list_of_x(), euler_solution.get_list_of_global_errors(self),
                  color=(255, 165, 0), name="euler GTE")
        self.plot(self.Max_Errors, euler_solution.get_n_list_of_error(),
                  euler_solution.get_list_of_max_global_errors(self), color=(255, 165, 0), name="euler GTE")

        # plot solution, GTE and Max GTE curves for ImprovedEulerMethod
        improved_euler_solution = ImprovedEulerMethod(self)
        self.plot(self.Solutions, improved_euler_solution.get_list_of_x(),
                  improved_euler_solution.get_list_of_y(improved_euler_solution.get_list_of_x()), color=(123, 104, 238),
                  name="improved euler")
        self.plot(self.Errors, improved_euler_solution.get_list_of_x(),
                  improved_euler_solution.get_list_of_global_errors(self), color=(123, 104, 238),
                  name="improved euler GTE")
        self.plot(self.Max_Errors, improved_euler_solution.get_n_list_of_error(),
                  improved_euler_solution.get_list_of_max_global_errors(self), color=(123, 104, 238),
                  name="improved euler GTE")

        # plot solution, GTE and Max GTE curves for RungeKuttaMethod
        runge_kutta_solution = RungeKuttaMethod(self)
        self.plot(self.Solutions, runge_kutta_solution.get_list_of_x(),
                  runge_kutta_solution.get_list_of_y(runge_kutta_solution.get_list_of_x()), color=(255, 255, 86),
                  name="runge-kutta")
        self.plot(self.Errors, runge_kutta_solution.get_list_of_x(),
                  runge_kutta_solution.get_list_of_global_errors(self), color=(255, 255, 86),
                  name="runge-kutta GTE")
        self.plot(self.Max_Errors, runge_kutta_solution.get_n_list_of_error(),
                  runge_kutta_solution.get_list_of_max_global_errors(self), color=(255, 255, 86),
                  name="runge-kutta GTE")

        # set checkboxes and radio button in default state
        self.GTE_button.setChecked(True)
        self.LTE_button.setChecked(False)
        self.exact_checkBox.setChecked(True)
        self.euler_checkBox.setChecked(True)
        self.improved_euler_checkBox.setChecked(True)
        self.runge_kutta_checkBox.setChecked(True)

    def clear_graphs(self):
        graphs = [self.Solutions, self.Errors, self.Max_Errors]
        for graph in graphs:
            graph.clear()
        app.processEvents()

    # to plot Solution graph for new input values
    def update_Solutions_plot(self, isExact, isEuler, isImproved_euler, isRunge_kutta):
        self.Solutions.clear()
        if isExact:
            exact_solution = ExactMethod(self)
            self.plot(self.Solutions, exact_solution.get_list_of_x(),
                      exact_solution.get_list_of_y(exact_solution.get_list_of_x()), color=(255, 20, 147), name="exact")
        if isEuler:
            euler_solution = EulerMethod(self)
            self.plot(self.Solutions, euler_solution.get_list_of_x(),
                      euler_solution.get_list_of_y(euler_solution.get_list_of_x()), color=(255, 165, 0), name="euler")
        if isImproved_euler:
            improved_euler_solution = ImprovedEulerMethod(self)
            self.plot(self.Solutions, improved_euler_solution.get_list_of_x(),
                      improved_euler_solution.get_list_of_y(improved_euler_solution.get_list_of_x()),
                      color=(123, 104, 238), name="improved euler")
        if isRunge_kutta:
            runge_kutta_solution = RungeKuttaMethod(self)
            self.plot(self.Solutions, runge_kutta_solution.get_list_of_x(),
                      runge_kutta_solution.get_list_of_y(runge_kutta_solution.get_list_of_x()), color=(255, 255, 86),
                      name="runge-kutta")
        app.processEvents()

    # for checkboxes
    def show_or_hide_curve(self):
        self.update_Solutions_plot(self.exact_checkBox.isChecked(), self.euler_checkBox.isChecked(),
                                   self.improved_euler_checkBox.isChecked(),
                                   self.runge_kutta_checkBox.isChecked())

    # to plot new graphs after choosing GTE or LTE curves
    def update_Errors_and_Max_Errors_plots(self, isShowGTE):
        self.Errors.clear()
        self.Max_Errors.clear()
        # if we choose GTE graphs
        if isShowGTE:
            euler_solution = EulerMethod(self)
            self.plot(self.Errors, euler_solution.get_list_of_x(), euler_solution.get_list_of_global_errors(self),
                      color=(255, 165, 0), name="euler GTE")
            self.plot(self.Max_Errors, euler_solution.get_n_list_of_error(),
                      euler_solution.get_list_of_max_global_errors(self), color=(255, 165, 0), name="euler LTE")
            improved_euler_solution = ImprovedEulerMethod(self)
            self.plot(self.Errors, improved_euler_solution.get_list_of_x(),
                      improved_euler_solution.get_list_of_global_errors(self), color=(123, 104, 238),
                      name="improved euler GTE")
            self.plot(self.Max_Errors, improved_euler_solution.get_n_list_of_error(),
                      improved_euler_solution.get_list_of_max_global_errors(self), color=(123, 104, 238),
                      name="improved euler GTE")
            runge_kutta_solution = RungeKuttaMethod(self)
            self.plot(self.Errors, runge_kutta_solution.get_list_of_x(),
                      runge_kutta_solution.get_list_of_global_errors(self), color=(255, 255, 86),
                      name="runge-kutta GTE")
            self.plot(self.Max_Errors, runge_kutta_solution.get_n_list_of_error(),
                      runge_kutta_solution.get_list_of_max_global_errors(self), color=(255, 255, 86),
                      name="runge-kutta GTE")
        # if we choose LTE graphs
        else:
            euler_solution = EulerMethod(self)
            self.plot(self.Errors, euler_solution.get_list_of_x(), euler_solution.get_list_of_local_errors(self),
                      color=(255, 165, 0), name="euler LTE")
            self.plot(self.Max_Errors, euler_solution.get_n_list_of_error(),
                      euler_solution.get_list_of_max_local_errors(self), color=(255, 165, 0), name="euler LTE")
            improved_euler_solution = ImprovedEulerMethod(self)
            self.plot(self.Errors, improved_euler_solution.get_list_of_x(),
                      improved_euler_solution.get_list_of_local_errors(self), color=(123, 104, 238),
                      name="improved euler LTE")
            self.plot(self.Max_Errors, improved_euler_solution.get_n_list_of_error(),
                      improved_euler_solution.get_list_of_max_local_errors(self), color=(123, 104, 238),
                      name="improved euler LTE")
            runge_kutta_solution = RungeKuttaMethod(self)
            self.plot(self.Errors, runge_kutta_solution.get_list_of_x(),
                      runge_kutta_solution.get_list_of_local_errors(self), color=(255, 255, 86),
                      name="runge-kutta LTE")
            self.plot(self.Max_Errors, runge_kutta_solution.get_n_list_of_error(),
                      runge_kutta_solution.get_list_of_max_local_errors(self), color=(255, 255, 86),
                      name="runge-kutta LTE")
        app.processEvents()

    # for radio buttons
    def show_or_hide_errors_curve(self):
        self.update_Errors_and_Max_Errors_plots(self.GTE_button.isChecked())


if __name__ == '__main__':
    main = MainWindow()
    main.show()
    app.processEvents()
    sys.exit(app.exec_())
