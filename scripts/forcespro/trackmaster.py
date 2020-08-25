"Based on example from https://codeloop.org/python-gui-how-to-create-paint-application-in-pyqt5/ by Parwiz"

from PyQt5.QtWidgets import QMainWindow, QApplication, QMenu, QMenuBar, QAction, QFileDialog
from PyQt5.QtGui import QIcon, QImage, QPainter, QPen, QBrush, QFont
from PyQt5.QtCore import Qt, QPoint
import sys
import InterpolateTrack
from python_sim_utils import plotter
import matplotlib.pyplot as plt
import numpy as np

class Window(QMainWindow):
    def __init__(self):
        super().__init__()


        title = "Paint Application"
        top = 400
        left = 400
        width = 800
        height = 600
        self.px_p_meter = 100

        icon = "icons/pain.png"

        self.setWindowTitle(title)
        self.setGeometry(top, left, width, height)
        self.setWindowIcon(QIcon(icon))

        self.image = QImage(self.size(), QImage.Format_RGB32)
        self.image.fill(Qt.white)

        self.waypoints = []
        self.drawing = False
        self.brushSize = 5
        self.brushColor = Qt.black
        self.lastPoint = QPoint()

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("File")
        brushSize = mainMenu.addMenu("Brush Size")
        brushColor = mainMenu.addMenu("Brush Color")
        '''
        saveAction = QAction(QIcon("icons/save.png"), "Save",self)
        saveAction.setShortcut("Ctrl+S")
        fileMenu.addAction(saveAction)
        saveAction.triggered.connect(self.save)
        '''
        exportAction = QAction(QIcon("icons/export.png"), "Export Track",self)
        exportAction.setShortcut("Ctrl+E")
        fileMenu.addAction(exportAction)
        exportAction.triggered.connect(self.exportTrack)

        clearAction = QAction(QIcon("icons/clear.png"), "Clear", self)
        clearAction.setShortcut("Ctrl+C")
        fileMenu.addAction(clearAction)
        clearAction.triggered.connect(self.clear)

        threepxAction = QAction( QIcon("icons/threepx.png"), "3px", self)
        brushSize.addAction(threepxAction)
        threepxAction.triggered.connect(self.threePixel)

        fivepxAction = QAction(QIcon("icons/fivepx.png"), "5px", self)
        brushSize.addAction(fivepxAction)
        fivepxAction.triggered.connect(self.fivePixel)

        sevenpxAction = QAction(QIcon("icons/sevenpx.png"),"7px", self)
        brushSize.addAction(sevenpxAction)
        sevenpxAction.triggered.connect(self.sevenPixel)

        ninepxAction = QAction(QIcon("icons/ninepx.png"), "9px", self)
        brushSize.addAction(ninepxAction)
        ninepxAction.triggered.connect(self.ninePixel)

        blackAction = QAction(QIcon("icons/black.png"), "Black", self)
        blackAction.setShortcut("Ctrl+B")
        brushColor.addAction(blackAction)
        blackAction.triggered.connect(self.blackColor)


        whitekAction = QAction(QIcon("icons/white.png"), "White", self)
        whitekAction.setShortcut("Ctrl+W")
        brushColor.addAction(whitekAction)
        whitekAction.triggered.connect(self.whiteColor)


        redAction = QAction(QIcon("icons/red.png"), "Red", self)
        redAction.setShortcut("Ctrl+R")
        brushColor.addAction(redAction)
        redAction.triggered.connect(self.redColor)

        greenAction = QAction(QIcon("icons/green.png"), "Green", self)
        greenAction.setShortcut("Ctrl+G")
        brushColor.addAction(greenAction)
        greenAction.triggered.connect(self.greenColor)

        yellowAction = QAction(QIcon("icons/yellow.png"), "Yellow", self)
        yellowAction.setShortcut("Ctrl+Y")
        brushColor.addAction(yellowAction)
        yellowAction.triggered.connect(self.yellowColor)

        painter = QPainter(self.image)
        painter.setBrush(Qt.black)
        painter.setPen(QPen(Qt.black, 2, Qt.SolidLine))
        painter.drawLine(10,50,10+self.px_p_meter,50)

        pen = QPen(Qt.black)
        pen.setWidth(2)
        painter.setPen(pen)

        font = QFont()
        font.setFamily('Times')
        font.setBold(True)
        font.setPointSize(12)
        painter.setFont(font)

        painter.drawText(20+self.px_p_meter, 55, "1 m")

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.drawing = True
            self.lastPoint = event.pos()
            self.waypoints.append([self.lastPoint.x(), self.lastPoint.y()])
            #print([self.lastPoint.x(), self.lastPoint.y()])
            #print(self.lastPoint)


    def mouseMoveEvent(self, event):
        if(event.buttons() & Qt.LeftButton) & self.drawing:
            painter = QPainter(self.image)
            painter.setPen(QPen(self.brushColor, self.brushSize, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin))
            painter.drawLine(self.lastPoint, event.pos())
            self.lastPoint = event.pos()
            self.update()



    def mouseReleaseEvent(self, event):

        if event.button() == Qt.LeftButton:
            self.drawing = False


    def paintEvent(self, event):
        canvasPainter  = QPainter(self)
        canvasPainter.drawImage(self.rect(),self.image, self.image.rect() )



    def save(self):
        filePath, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG(*.png);;JPEG(*.jpg *.jpeg);;All Files(*.*) ")
        fig, ax1= plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))
        points = 1/100 * np.array(self.waypoints)

        ax1.scatter(points[:,0], points[:,1])
        if filePath == "":
            return
        self.image.save(filePath)

    def exportTrack(self):
        print("exporting Track")
        fig, ax1= plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))
        points = 1/self.px_p_meter*np.array(self.waypoints)
        points[:,0] = points[:,0] - np.mean(points[:,0])
        points[:,1] = points[:,1] - np.mean(points[:,1])
        ax1.scatter(points[:,0], points[:,1])
        plt.show()
        name = QFileDialog.getSaveFileName(self, "Save Waypoints", "", "CSV(*.csv)")
        np.savetxt(name[0], points, delimiter = ', ')
        track_lu_table, smax = InterpolateTrack.generatelookuptable("tracks/"+name[0][53:-4])
        r = 0.2
        lencar = 0.06
        trk_plt = plotter(track_lu_table, smax, r, lencar)
        trk_plt.plot_track()
        plt.show()

    def clear(self):
        print("clearing")
        self.image.fill(Qt.white)
        self.waypoints = []
        self.update()


    def threePixel(self):
        self.brushSize = 3

    def fivePixel(self):
        self.brushSize = 5

    def sevenPixel(self):
        self.brushSize = 7

    def ninePixel(self):
        self.brushSize = 9


    def blackColor(self):
        self.brushColor = Qt.black

    def whiteColor(self):
        self.brushColor = Qt.white

    def redColor(self):
        self.brushColor = Qt.red

    def greenColor(self):
        self.brushColor = Qt.green

    def yellowColor(self):
        self.brushColor = Qt.yellow




if __name__ == "__main__":
    print("=============== Welcome to Trackmaster ===============")
    print("This is a GUI for drawing tracks to use with the MPCC controller.")
    print("Use this application to draw a sparse waypoint trajectory and hit Ctrl+ E to generate the dense lookup table with arclength parametrization.")
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    app.exec()
