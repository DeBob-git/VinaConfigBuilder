import sys 
import os

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QDesktopWidget, QLabel, QPushButton, QInputDialog, QLineEdit, QFileDialog, QAction
from PyQt5.QtWidgets import QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget

from PyQt5.QtGui import QIcon
from PyQt5 import Qt

from PyQt5 import QtGui
from PyQt5 import QtCore

try:
    import plip
    from plip.structure.preparation import PDBComplex
except:
    print("No plip in the system")

import frst2  


# ==========================================================================        

 
class atm:
    def __init__(self):
        self.tip=''
        self.letter=''
        self.pos = None
        self.lname=''
        self.sname=''
        self.X = None
        self.Y = None
        self.Z = None
        self.selected=False
        self.amks=[]

class point:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.X=x
        self.Y=y
        self.Z=z
      
class cube:
    def __init__(self, cp=point(0,0,0), xl=1 ,yl=1,zl=1):
        self.CP = cp
        self.xl=xl
        self.yl=yl
        self.zl=zl
        self.shape = []       
        
    def _fill_shapes(self):
        for i in range(0,8):
            self.shape.append(point())         

        self.shape[0]=[(self.CP.X-self.xl/2),(self.CP.Y-self.yl/2),(self.CP.Z-self.zl/2)] # 1
        self.shape[1]=[(self.CP.X-self.xl/2),(self.CP.Y+self.yl/2),(self.CP.Z-self.zl/2)] # 2
        self.shape[2]=[(self.CP.X+self.xl/2),(self.CP.Y-self.yl/2),(self.CP.Z-self.zl/2)] # 3
        self.shape[3]=[(self.CP.X-self.xl/2),(self.CP.Y-self.yl/2),(self.CP.Z+self.zl/2)] # 4
        
        self.shape[4]=[(self.CP.X+self.xl/2),(self.CP.Y+self.yl/2),(self.CP.Z-self.zl/2)] # 5
        self.shape[5]=[(self.CP.X+self.xl/2),(self.CP.Y-self.yl/2),(self.CP.Z+self.zl/2)] # 6  
        self.shape[6]=[(self.CP.X+self.xl/2),(self.CP.Y+self.yl/2),(self.CP.Z+self.zl/2)]
        self.shape[7]=[(self.CP.X-self.xl/2),(self.CP.Y-self.yl/2),(self.CP.Z+self.zl/2)]      
      
# ======================================================================================================================

class mainApp(QtWidgets.QMainWindow,frst2.Ui_MainWindow):
    def __init__(self):

        super().__init__()
        self.setupUi(self)  
               
        
        self.exitButt.clicked.connect(self._exit)
        self.LIGC_BUTT.clicked.connect(self._ligCentroid)
        self.LIGC_BUTT_PLIP.clicked.connect(self._ligPlip)
        
        self.AMKC_BUTT.clicked.connect(self._amkCentroid)
        
        #self.autoCalc_CNT.clicked.connect(lambda checked,tip=2: self._amkCentroid(tip))
        self.autoCalc_AMK.clicked.connect(lambda checked,tip=1: self._autoCalc(tip))
        self.autoCalc_LIG.clicked.connect(lambda checked,tip=2: self._autoCalc(tip))
        self.autoCalc.clicked.connect(lambda checked,tip=4: self._autoCalc(tip))
        self.autoCalc_PL_L.clicked.connect(lambda checked,tip=3: self._autoCalc(tip))
        
        self.allIn.itemClicked.connect(self._allInChanged)        
        
        self.getFileButt.clicked.connect(self._fileSelect)
        self.getRecFileButt.clicked.connect(self._recFileSelect)
        
        self.saveButt.clicked.connect(self._saveResults)
        
        self.sInAMK.textChanged.connect(self._calcAMK)
        
        self.stepSides.valueChanged.connect(self._sidesStep)
        self.stepCenter.valueChanged.connect(self._centerStep)
        
        self.xl.valueChanged.connect(self._calcAMK)
        self.yl.valueChanged.connect(self._calcAMK)
        self.zl.valueChanged.connect(self._calcAMK)
        
        self.xC.valueChanged.connect(self._calcAMK)
        self.yC.valueChanged.connect(self._calcAMK)
        self.zC.valueChanged.connect(self._calcAMK)  
        
       # self.allIn.itemSelectionChanged.connect(self._allInSelectionChanged)
            
    
        self.atms=[]
        self.atmDict={}
        self.myCube=cube(point(0,0,0),50,50,50)
        self.myCube._fill_shapes()
        self.currLIG.setText('None')
        self.gamma.setText('0')
        
# ==================================================================================================

    def _allInChanged(self):
        print('Pressed ' + str(self.allIn.currentRow()) + self.allIn.currentItem().text() )
        if self.allIn.currentRow()!=-1:
            key=self.allIn.currentItem().text()
            kkk=key.split('_')
            if len(kkk)==1:
                print('Not a ligand')
                self.currLIG.setText('None')
            else:
                self.currLIG.setText(key)

# ==================================================================================================
   
    def _saveResults(self):
        if self.fileToScan.text()=='' or self.fileToScan.text()==None:
            shortName='save.txt'
            sName='nopdb.pdb'
        else:
            shortName=os.path.splitext(self.fileToScan.text())[0] + '.txt'
            sName=os.path.basename(self.fileToScan.text())
            
        fileName=QFileDialog.getSaveFileName(self,"Save CUBE params",shortName,'*.txt')
        print(fileName[0])
        
        if fileName[0]!='':
            f = open(fileName[0], "w")
            #f.write('receptor=' + self.recept.text())
            #f.write('\n')
            f.write('exhaustiveness=' + str(self.exhaust.value()))
            f.write('\n')
            
            f.write('center_x=' + str(self.xC.value()) )
            f.write('\n')
            f.write('center_y=' + str(self.yC.value()) )
            f.write('\n')
            f.write('center_z=' + str(self.zC.value()) )
            f.write('\n')
            
            f.write('size_x=' + str(self.xl.value()))
            f.write('\n')
            f.write('size_y=' + str(self.yl.value()))
            f.write('\n')
            f.write('size_z=' + str(self.zl.value()))
            f.write('\n')
            
            f.write('num_modes=' + str(self.num_modes.value()))
            f.write('\n')            
            
            f.close()
 
 # ==================================================================================================                  

    def _sidesStep(self):
        self.xl.setSingleStep(self.stepSides.value())
        self.yl.setSingleStep(self.stepSides.value())
        self.zl.setSingleStep(self.stepSides.value())
  
    def _centerStep(self):
        self.xC.setSingleStep(self.stepCenter.value())
        self.yC.setSingleStep(self.stepCenter.value())
        self.zC.setSingleStep(self.stepCenter.value())        

# ==================================================================================================
       
    def _showItem(self):
        ss=self.plipAll.Item()
        print (ss.Text)

# ==================================================================================================
        
    def _getCubeParams(self):
        xcc=self.xC.value()
        ycc=self.yC.value()
        zcc=self.zC.value()      
        #print(xcc,ycc,zcc)
        
        xl=self.xl.value()
        yl=self.yl.value()
        zl=self.zl.value()        
        #print(xl,yl,zl)
        
        self.myCube.CP=point(xcc,ycc,zcc)
        self.myCube.xl=xl
        self.myCube.yl=yl
        self.myCube.zl=zl
        
        self.myCube._fill_shapes()

# ==================================================================================================        
          
    def _exit(self):
        print("See your later.")
        self.close()

# ==================================================================================================
        
    def _fileSelect(self):
        self.fileToScan.setText(QFileDialog.getOpenFileName(self,'Select PDB file',None,'*.pdb')[0])
        if self.fileToScan.text()!='':
            self._clear_list()
            self.sInAMK.setText('')
            self._fill_AMK_LW()
            self._calcAMK()

# ==================================================================================================
            
    def _recFileSelect(self):
        recFileDir=QFileDialog.getOpenFileName(self,'Select receptor file (pdbqt)',None,'*.pdbqt')[0]
        if recFileDir!='':
            recFile=os.path.basename(recFileDir)
            self.recept.setText(recFile)

# ==================================================================================================        
    
    def _clear_list(self):
        self.partAMK.clear()
        self.fullAMK.clear()                
        self.partLIG.clear()
        self.fullLIG.clear()   
        self.allIn.clear()

# ==================================================================================================
        
    def _clear_coord(self):     
        self.xC.setValue(0)
        self.yC.setValue(0)
        self.zC.setValue(0)        
        
        self.xl.setValue(50)
        self.yl.setValue(50)
        self.zl.setValue(50)                

# ==================================================================================================
        
    def _calcAMK(self):
        
        self._getCubeParams()
        self._clear_list()
        
        sPlip = self.sInAMK.text()
        sPlip = sPlip.split(',')
        errs=0
# ============================= fill AMK        
        for sp in sPlip:
            if not sp in self.atmDict:
                #print('No such AMK in PDB ' + sp)
                errs+=1
                #if sp != '' and sp != None:
                #    buttonReply = QMessageBox.question(self, 'Warning!!', "No such AMK " + sp ,QMessageBox.Ok)
            else:
                atmElement=self.atmDict[sp]
                line = sp
                sname = atmElement[0].sname
                tip=atmElement[0].tip
                
                inCube = self._inCube(atmElement)
                item=Qt.QListWidgetItem(line)
                self.allIn.addItem(item) 
                
                item.setBackground(QtGui.QColor('pink'))

                if inCube==1:
                    if tip=='ATOM':
                        self.partAMK.addItem(line)  
                        item.setBackground(QtGui.QColor('yellow'))
                if inCube==2:
                    if tip=='ATOM':
                        self.fullAMK.addItem(line)
                        item.setBackground(QtGui.QColor('lime'))
                        
                        
# ============================= fill LIG                                
        for key in self.atmDict:
            atmElement=self.atmDict[key]
            sname = atmElement[0].sname
            lname = atmElement[0].lname
            tip=atmElement[0].tip
                        
            if atmElement[0].tip=='HETATM':
                item=Qt.QListWidgetItem(sname+'_'+lname)
                self.allIn.addItem(item) 
                item.setBackground(QtGui.QColor('pink'))
                inCube = self._inCube(atmElement)
                if inCube==1:
                    self.partLIG.addItem(sname)   
                    item.setBackground(QtGui.QColor('yellow'))
                if inCube==2:
                    self.fullLIG.addItem(sname)  
                    item.setBackground(QtGui.QColor('lime'))

                #self.allIn.addItem(sname)        
                
# ==================================================================================================                           
                    
    def _incube_check(self,cube,point):
        bresx = (float(cube.CP.X)-float(cube.xl/2)) <= float(point.X) and (float(point.X) <= (float(cube.CP.X)+float(cube.xl/2)))
        bresy = (float(cube.CP.Y)-float(cube.yl/2)) <= float(point.Y) and (float(point.Y) <= (float(cube.CP.Y)+float(cube.yl/2)))
        bresz = (float(cube.CP.Z)-float(cube.zl/2)) <= float(point.Z) and (float(point.Z) <= (float(cube.CP.Z)+float(cube.zl/2)))
       #print (bresx, bresy,bresz)
        return (bresx and bresy and bresz)
    
# ==================================================================================================

    def _inCube(self,AMK):
       # print(len(AMK))
       # self._getCubeParams()
        
        inCube=0
        for x in AMK:
           #for t in am:
           #print(x.sname,x.pos,x.X,x.Y,x.Z)
           
           if self._incube_check(self.myCube,point(x.X,x.Y,x.Z)):
               inCube+=1
               
        if inCube==len(AMK):
            return 2
        if inCube>0 and inCube<len(AMK):
            return 1
        if inCube==0:
            return 0                

# ==================================================================================================

    def _ligCentroid(self):
        X=[]
        Y=[]
        Z=[]
        
        key=self.currLIG.text()
        if key=='None':
            return 0,0,0
        
        kkk=key.split('_')
        if len(kkk)==1:
            print('Not a ligand')
            return -1
         
        atmElement=self.atmDict[kkk[1]]
        for x in atmElement:
            X.append(float(x.X))
            Y.append(float(x.Y))
            Z.append(float(x.Z))
         
        centr=sum(X)/len(X),sum(Y)/len(Y),sum(Z)/len(Z)
        self.myCube.CP.X=sum(X)/len(X)
        self.myCube.CP.y=sum(Y)/len(Y)
        self.myCube.CP.z=sum(Z)/len(Z)
        
        self.xC.setValue(sum(X)/len(X))
        self.yC.setValue(sum(Y)/len(Y))
        self.zC.setValue(sum(Z)/len(Z))

        print(centr)
        return centr    
    
# ==================================================================================================    

    def _amkCentroid(self,tip=1):
        X=[]
        Y=[]
        Z=[]
        
        sPlip = self.sInAMK.text()
        sPlip = sPlip.split(',')
        
        for sp in sPlip:
            if sp in self.atmDict:
                atmElement=self.atmDict[sp]  
                for x in atmElement:
                    X.append(float(x.X))
                    Y.append(float(x.Y))
                    Z.append(float(x.Z))               
                    
        if tip==2:
            for key in self.atmDict:
                atmElement=self.atmDict[key]
                if atmElement[0].tip=='HETATM':
                    for x in atmElement:
                        X.append(float(x.X))
                        Y.append(float(x.Y))
                        Z.append(float(x.Z))  
                        
        if tip==3:
            for key in self.atmDict:
                atmElement=self.atmDict[key]
                
                for x in atmElement:
                    X.append(float(x.X))
                    Y.append(float(x.Y))
                    Z.append(float(x.Z))        
        
        if len(X)!=0 and len(Y)!=0 and len(Z)!=0:                
            self.xC.setValue(sum(X)/len(X))
            self.yC.setValue(sum(Y)/len(Y))
            self.zC.setValue(sum(Z)/len(Z))         
        
# ====================================================================================    
    
    def _autoCalc(self,tip):
        X=[]
        Y=[]
        Z=[]  
        
        
        xC=self.xC.value()
        yC=self.yC.value()
        zC=self.zC.value()
        gamma=float(self.gamma.text())
        
        
        # ================================ all AMK from PLIP string
        sPlip = self.sInAMK.text()
        sPlip = sPlip.split(',')
        
        if tip==1 or tip==3: 
            for sp in sPlip:
                if sp in self.atmDict:
                    atmElement=self.atmDict[sp]  
                    for x in atmElement:
                        X.append(abs((float(x.X)-xC)))
                        Y.append(abs((float(x.Y)-yC)))
                        Z.append(abs((float(x.Z)-zC)))
                        
        # ================================ Just selected LIG
        if tip==2:
            key=self.currLIG.text()
            if key!='None':
                kkk=key.split('_')
                if len(kkk)==1:
                    print('Not a ligand')
                    return -1
                atmElement=self.atmDict[kkk[1]]                
                if atmElement[0].tip=='HETATM':
                    for x in atmElement:
                        X.append(abs((float(x.X)-xC)))
                        Y.append(abs((float(x.Y)-yC)))
                        Z.append(abs((float(x.Z)-zC)))   
                        
        # ================================ all AMK and selected LIG  
        
        if tip==4:
            self._amkCentroid(3)
            for key in self.atmDict:
                atmElement=self.atmDict[key]
                for x in atmElement:
                    X.append(abs((float(x.X)-xC)))
                    Y.append(abs((float(x.Y)-yC)))
                    Z.append(abs((float(x.Z)-zC)))                         
                        
        # ================================ all AMK and all selected LIG                
        if tip==3:
            for key in self.atmDict:
                atmElement=self.atmDict[key]
                if atmElement[0].tip=='HETATM':
                    for x in atmElement:
                        X.append(abs((float(x.X)-xC)))
                        Y.append(abs((float(x.Y)-yC)))
                        Z.append(abs((float(x.Z)-zC)))     
                        
                    
        if len(X)>0 and len(Y)>0 and len(Z)>0:
            self.xl.setValue(round(max(X)*2,2)+0.01+gamma)
            self.yl.setValue(round(max(Y)*2,2)+0.01+gamma)
            self.zl.setValue(round(max(Z)*2,2)+0.01+gamma)
            return max(X),max(Y),max(Z)           

# ==================================================================================================
    
    def _fill_AMK_LW(self):
        
        try:
            my_mol = PDBComplex()
            my_mol.load_pdb(self.fileToScan.text())
            my_mol.analyze()
        except:
            print("Can't load plip functions.")
               
        
        self._clear_list()
        self._clear_coord()
        self.currLIG.setText('None')
        self.atmDict={}
        amkLig=[]
        
        f=open(self.fileToScan.text(),'r')
        amkCnt=0
        ligCnt=0
        for st in f:
            lst=st.split(' ')
            lst =' '.join(st.split())
            lst=lst.split(' ')
            
            if lst[0][0:4] == 'ATOM' or lst[0][0:6] == 'HETATM':
                tkatm = atm()
                if lst[0][0:4] == 'ATOM':
                    tkatm.tip = lst[0][0:4]
                    
                if len(lst[0])>=6 and lst[0][0:6] == 'HETATM':
                    lst_n=st[:6] + ' ' + st[6:]
                    lst=' '.join(lst_n.split())
                    lst=lst.split(' ')
                    tkatm.tip = lst[0]
                    
                tkatm.pos=lst[1]
                tkatm.sname=lst[3]
                tkatm.letter = lst[4]
                
                if len(lst[4])>1:
                    tkatm.lname=lst[4]
                    tkatm.X=lst[5]
                    tkatm.Y=lst[6]
                    tkatm.Z=lst[7]                    
                else:
                    tkatm.lname=lst[4]+lst[5]
                    tkatm.X=lst[6]
                    tkatm.Y=lst[7]
                    tkatm.Z=lst[8]                    

                if tkatm.lname in self.atmDict:
                    xx = []
                    xx = self.atmDict[tkatm.lname]
                    xx.append(tkatm)
                    self.atmDict[tkatm.lname] = xx
                    #print(xx)
                    # ================='Just add one more atom'
                else:
                    xx = []
                    xx.append(tkatm)                    
                    self.atmDict[tkatm.lname]=xx
                    #print(xx)
                    # ================='New item into DICT atm'
                    
        allAmk=0            
        for key in self.atmDict:
            atmElement=self.atmDict[key]
            seen=set()
            
            if atmElement[0].tip == 'HETATM': 
                amklig=[]
                my_bsid=atmElement[0].sname+':'+atmElement[0].lname[0]+':'+atmElement[0].lname[1:]
                try:
                    my_interactions = my_mol.interaction_sets[my_bsid] 
                             
                    print(my_bsid) 
                
                    #for cnts in my_interactions.all_hydrophobic_contacts:
                    for cnts in my_interactions.interacting_res:
                        #alig=cnts[7]+str(cnts[6])
                        alig=cnts[-1:]+cnts[:-1]
                        if not (alig in seen):
                            amklig.append(alig)
                            allAmk+=1
                            seen.add(alig)
                    atmElement[0].amkLig=amklig    
                    print (amklig) 
                except:
                    print('No such lig in pdb')
                    atmElement[0].amkLig=amklig 
                    
        print('all AMK :' + str(allAmk))        

# ==================================================================================================
                
    def _ligPlip(self):
        key=self.currLIG.text()
        if key!='None':
            kkk=key.split('_')
            if len(kkk)==1:
               print('Not a ligand')
               return -1
            atmElement=self.atmDict[kkk[1]]
            print(atmElement[0].amkLig)
            self.sInAMK.setText(','.join(atmElement[0].amkLig)  )
                
# A99,A101,A201,A206,A208,A245,A64,A203,A204,A215,A216
# =====================================================================================  

def main():
    app = QtWidgets.QApplication(sys.argv)  # РќРѕРІС‹Р№ СЌРєР·РµРјРїР»СЏСЂ QApplication
    window = mainApp()                      # РЎРѕР·РґР°С‘Рј РѕР±СЉРµРєС‚ РєР»Р°СЃСЃР° ExampleApp
    window.show()                           # РџРѕРєР°Р·С‹РІР°РµРј РѕРєРЅРѕ
    app.exec_()                             # Рё Р·Р°РїСѓСЃРєР°РµРј РїСЂРёР»РѕР¶РµРЅРёРµ

if __name__ == '__main__':  
    main()  