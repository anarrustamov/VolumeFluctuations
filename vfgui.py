#Author: A.Rustamov 

from importlib.metadata import distribution
from turtle import st
from volfluct import Volfluct
from vfguibase import Ui_Form
from qtpy.QtWidgets import QTextEdit
import sympy as sp

class VF(Ui_Form):
    def __init__(self, Form):
        self.setupUi(Form)
        self.pushButton.clicked.connect(self.Derive)
        self.comboBox_DeriveFor.addItem("Volume fluctuations")
        self.comboBox_DeriveFor.addItem("Correction for Volume fluctuations")
    def Derive(self):
        cumorder = self.spinBox_cumorder.value()
        a = Volfluct()
        results = []
        names = []
        mixresults = []
        mixnames = []
        value = str(self.comboBox_DeriveFor.currentText())
        if value == "Volume fluctuations":
            self.clearTextBrowsers()
            for i in range(cumorder):
                [r, n] = a.getCumulant1(i+1, '', ['n','N'])
                results.extend([r])
                names.extend([n])
            self.textBrowser.append("Deriving formulas for volume fluctuations for order: " + str(cumorder))
            self.textBrowser.append("__________________________________________")
           
            self.AddDescription1()
            for i in range(len(results)):
                self.textBrowser.append(str(names[i]) + " = " + str(results[i]))
                self.textBrowser.append(" ")

            if self.radioButton.isChecked():
                for i in range(cumorder):
                    for j in range(cumorder):
                        if (i+j+2) > cumorder:
                            continue
                        if ( i > j):
                            continue
                        [r, n] = a.getCrossCumulant1(i+1,j+1, '', ['n1','N1'], ['n2', 'N2'])
                        mixresults.extend([r])
                        mixnames.extend([n])
                self.textBrowser_2.append("Deriving formulas for volume fluctuations for order: " + str(cumorder))
                self.textBrowser_2.append("__________________________________________")
                for i in range(len(mixresults)):
                    self.textBrowser_2.append(str(mixnames[i]) + " = " + str(mixresults[i]))
                    self.textBrowser_2.append(" ")
                self.AddDescription2()
        if value == "Correction for Volume fluctuations":
            self.clearTextBrowsers()
            self.textBrowser.append("Deriving formulas for correcting volume fluctuations for order: " + str(cumorder))
            self.textBrowser.append("__________________________________________")
            if self.radioButton.isChecked():
                [results1, resultsname1] = a.getCumW2(cumorder, 1)
            else:
                [results1, resultsname1] = a.getCumW2(cumorder, 0)
            self.textBrowser.append("corrected cumulants:")
            self.textBrowser.append("")
            for i in range(1, len(results1[0])):
                self.textBrowser.append(str(resultsname1[0][i]) + " = " + str(results1[0][i]))
                self.textBrowser.append(" ")
            self.textBrowser.append("bias terms")
            self.textBrowser.append("")
            for i in range(1, len(resultsname1[1])):
                self.textBrowser.append(str(resultsname1[1][i]) + " =  " + str(results1[1][i]))
                self.textBrowser.append(" ")
            self.AddDescription3()
            if self.radioButton.isChecked():
                self.textBrowser_2.setAcceptRichText(True)
                self.textBrowser_2.append("Deriving formulas for correcting volume fluctuations for order: " + str(cumorder))
                self.textBrowser_2.append("__________________________________________")
                self.textBrowser_2.append("corrected cumulants:")
                self.textBrowser_2.append("")
                for i in range(0, len(results1[2])):
                    self.textBrowser_2.append(str(resultsname1[2][i]) + " = " + str(results1[2][i]))
                    self.textBrowser_2.append(" ")
                self.textBrowser_2.append("__________________________________________")
                self.textBrowser_2.append("bias terms")
                self.textBrowser_2.append(" ")

                for i in range(0, len(resultsname1[3])):
                    self.textBrowser_2.append(str(resultsname1[3][i]) + " = " + str(results1[3][i]))
                    self.textBrowser_2.append(" ")
                
                self.AddDescription4()

    def clearTextBrowsers(self):
        self.textBrowser.clear()
        self.textBrowser.clearHistory()
        self.textBrowser_3.clear()
        self.textBrowser_3.clearHistory()
        self.textBrowser_2.clear()
        self.textBrowser_2.clearHistory()

    def AddDescription1(self):
        self.textBrowser_3.clear()
        self.textBrowser_3.clearHistory()
        self.textBrowser_3.append("Description: Pure cumulants, including volume fluctuations")
        self.textBrowser_3.append("__________________________________________")
        self.textBrowser_3.append(" ")
        self.textBrowser_3.append("k_n[N] -> nth order cumulant of the multiplicity distribution (can be measured)")
        self.textBrowser_3.append("k_n[n] -> nth order cumulant of the multiplicity distribution per single source") 
        self.textBrowser_3.append("<W> -> mean number of sources") 
        self.textBrowser_3.append("<n> -> mean number of particles per source") 

    def AddDescription2(self):
        self.AddDescription1()
        self.textBrowser_3.append(" ")
        self.textBrowser_3.append("Description: Mixed cumulants, including volume fluctuations")
        self.textBrowser_3.append("__________________________________________")
        self.textBrowser_3.append(" ")
        self.textBrowser_3.append("k_mn[N1N2] -> covariance of order mn between multiplicity distributions of particles N1 and N2")
        self.textBrowser_3.append("k_mn[n1n2] -> covariance of order mn between multiplicity distributions of particles n1 and n2 (per source)")
        
    def AddDescription3(self):
        self.textBrowser_3.clear()
        self.textBrowser_3.clearHistory()
        self.textBrowser_3.append("Description: Correction formulas for pure cumulants")
        self.textBrowser_3.append("__________________________________________")
        self.textBrowser_3.append(" ")
        self.textBrowser_3.append("k_corr[N] -> corrected nth order cumulant against volume fluctuations")
        self.textBrowser_3.append("k_n[N] -> nth order cumulant of multiplicity distribution (measured)")
        self.textBrowser_3.append("C_n[M] -> nth order factorial cumulant of total multiplicity distribution (measured)")
        self.textBrowser_3.append("delta[n] -> bias term for the nth order pure cumulant")
        self.textBrowser_3.append("C_bar_n[M] -> nth order factorial cumulant of total multiplicity distribution for a fixed volume")
        
    def AddDescription4(self):
        self.AddDescription3()
        self.textBrowser_3.append(" ")
        self.textBrowser_3.append("Description: Correction formulas for mixed cumulants")
        self.textBrowser_3.append("__________________________________________")
        self.textBrowser_3.append(" ")
        self.textBrowser_3.append("k_nm_corr[N1N2] -> corrected covariance of order mn between multiplicity distributions of particles N1 and N2")
        self.textBrowser_3.append("k_nm_[N1N2] -> covariance of order mn between multiplicity distributions of particles N1 and N2 (measured)")
        self.textBrowser_3.append("Delta[nm] -> bias term covariance of order mn between multiplicity distributions of particles N1 and N2")





