from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt


class EbeamTable():

    def __init__(self):
        self.init_table()


    def init_table(self):

        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(2)
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        self.table.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)


    def update_table(self, ebp):

        if ebp[0] is None:
            self.table.setRowCount(1)
            self.table.setSpan(0, 0, 1, 2)

            item = QtWidgets.QTableWidgetItem("No solution")
            item.setTextAlignment(Qt.AlignCenter)
            self.table.setItem(0, 0, item)

            return

        n = len(ebp) - 1
        self.table.setRowCount(n)
        self.table.clearSpans()
        
        for i in range(n):
            item = QtWidgets.QTableWidgetItem(ebp[i+1][0])
            self.table.setItem(i, 0, item)

            item = QtWidgets.QTableWidgetItem(str(ebp[i+1][1]))
            item.setTextAlignment(Qt.AlignCenter)
            self.table.setItem(i, 1, item)
