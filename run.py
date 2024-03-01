from vfgui import VF
from qtpy import QtWidgets
import sys

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = VF(Form)
    Form.show()
    sys.exit(app.exec())

