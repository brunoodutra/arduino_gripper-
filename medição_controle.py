from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
from pyqtgraph.ptime import time
import serial

app = QtGui.QApplication([])

win = pg.GraphicsWindow()
win.setWindowTitle('pyqtgraph example: logAxis')

# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)

p1 = win.addPlot(1,1, title="u_GMV - sinal de controle")
p2 = win.addPlot(1,0, title="Y - sinal medido e referência")
#p3 = win.addPlot(2,0, title="Yr - referencia")
p4 = win.addPlot(1,2, title="e - erro")

p5 = win.addPlot(0,0, title="flexor ")
p6 = win.addPlot(0,1, title="extensor")
p7 = win.addPlot(0,2, title="movimento estimado")


p1.showGrid(True, True)
p2.showGrid(True, True)
#p3.showGrid(True, True)
p4.showGrid(True, True)
p5.showGrid(True, True)
p6.showGrid(True, True)
p7.showGrid(True, True)

p1.setRange(xRange=[0,200])
p2.setRange(xRange=[0,200])
p4.setRange(xRange=[0,200])

p5.setRange(xRange=[0,1000])
p6.setRange(xRange=[0,1000])
p7.setRange(xRange=[0,1000])

#p1.setDownsampling(mode='peak')
#p2.setDownsampling(mode='peak')
#p4.setDownsampling(mode='peak')

#p3.setClipToView(True)
#p4.setClipToView(True)
#p3.setRange(xRange=[1000,0])
#p3.setLimits(xMax=1000)

#p1.setLabel('bottom', 'Index', units='B')

curve1 = p1.plot()
curve2 = p2.plot(name="curve1")
curve3 = p2.plot(name="curve2")
curve4 = p4.plot()

curve5 = p5.plot(name="flexor1")
curve5_1 = p5.plot(name="fk1")
curve6 = p5.plot(name="flexor2")
curve6_1 = p5.plot(name="fk2")
curve7 = p5.plot(name="fusao1")

curve8 = p6.plot(name="extensor1")
curve8_1 = p6.plot(name="fk3")
curve9 = p6.plot(name="extensor2")
curve9_1 = p6.plot(name="fk4")
curve10 = p6.plot(name="fusao2")

curve11 = p7.plot()


data1 = [0];
data2 = [0];
data3 = [0];
data4 = [0];
data5 = [0];
data5_1 = [0];
data6 = [0];
data6_1 = [0];
data7 = [0];
data8 = [0];
data8_1 = [0];
data9 = [0];
data9_1 = [0];
data8 = [0];
data9 = [0];
data10 = [0];

data11 = [0];

raw=serial.Serial('COM4',250000)
arquivo2 = open('movimentos.txt','w')
arquivo = open('controle.txt','w')

ptr = 0; ptr2 = 0
d=0; d2=0;
ok=0; ok2=0;
init=0;
lastTime = time()
fps = None
line="  ";
def Serial():
    global line;
    line = raw.readline()
def update():
    
    global line,init,ok,curve1,cuve2,curve3, data,ptr,data1,data2,data3,data4,d,lastTime, fps
    #p1.setTitle('Emg data 1   %0.2f Hz' % fps)
    #print(line)   
    if ("T" in line):
       init=1;
       arquivo.write(str(line[0:len(line)-4]+'\r\n'));
       line=line.split(",")
       if len(line)>3:
            data1.append(float(line[0]))
            data2.append(float(line[1]))
            data3.append(float(line[2]))
            data4.append(float(line[3]))
            ok=1;
            app.processEvents()
            ptr=ptr+1;
            #app.processEvents()  ## force complete redraw for every plot
    else:
        if(init==0):
            raw.write("s");#GMV
            #raw.write("d");#PID
            1;
    if(ptr>199):                            #If you have 50 or more points, delete the first one from the array
            #xdata1=xdata1*0;
            data1=data1*0;
            data2=data2*0;
            data3=data3*0;
            data4=data4*0;
            ptr=0;
            d=d+200;
            #p2.setXRange(d, d+500, padding=0)
    
def update2():
    now=time()
    global lastTime,line,init,ok2,curve5,cuve6,curve7, data,ptr2,data5,data6,data5_1,data6_1,data7,data8,data9,data8_1,data9_1,data10,data11,d2
    dt = now - lastTime
    lastTime = now        
    fps = 1/dt
    #print(fps)
    #p1.setTitle('Emg data 1   %0.2f Hz' % fps)
    if ("E" in line):
        
       arquivo2.write(str(line[0:len(line)-4]+'\r\n'));
       line=line.split(",")
       if len(line)>10:
            data5.append(float(line[0]))
            data5_1.append(float(line[1]))
            data6.append(float(line[2]))
            data6_1.append(float(line[3]))
            data7.append(float(line[4]))
            data8.append(float(line[5]))
            data8_1.append(float(line[6]))
            data9.append(float(line[7]))
            data9_1.append(float(line[8]))
            data10.append(float(line[9]))
            data11.append(float(line[10]))
            ok2=1;
            app.processEvents()
            ptr2=ptr2+1;

    if(ptr2>999):                            #If you have 50 or more points, delete the first one from the array
            #xdata1=xdata1*0;
            data5=data5*0;
            data5_1=data5_1*0;
            data6=data6*0;
            data6_1=data6_1*0;
            data7=data7*0;
            data8=data8*0;
            data8_1=data8_1*0;
            data9=data9*0;
            data9_1=data9_1*0;
            data10=data10*0;
            data11=data11*0;
            ptr2=0;
            d2=d2+1000;
            #p2.setXRange(d, d+500, padding=0)    
def plots():
    global ok,ok2,data1,xdata2,xdata3,xdata4,data5,xdata6,xdata7,data8,xdata9,xdata10,xdata11
    
    if(ok==1):
                xdata1 = np.array(data1,dtype='float64')
                xdata2 = np.array(data2,dtype='float64')
                xdata3 = np.array(data3,dtype='float64')
                xdata4 = np.array(data4,dtype='float64')
                curve1.setData(data1,pen='w')
                curve2.setData(xdata2,pen='b')
                curve3.setData(xdata3,pen='r')
                curve4.setData(xdata4,pen='y')
                ok=0;

    if(ok2==1):
                xdata5 = np.array(data5,dtype='float64')
                xdata5_1 = np.array(data5_1,dtype='float64')
                xdata6 = np.array(data6,dtype='float64')
                xdata6_1 = np.array(data6_1,dtype='float64')
                xdata7 = np.array(data7,dtype='float64')
                xdata8 = np.array(data8,dtype='float64')
                xdata8_1 = np.array(data8_1,dtype='float64')
                xdata9 = np.array(data9,dtype='float64')
                xdata9_1 = np.array(data9_1,dtype='float64')
                xdata10 = np.array(data10,dtype='float64')
                xdata11 = np.array(data11,dtype='float64')
                curve5.setData(xdata5,pen='w')
                curve5_1.setData(xdata5_1,pen='m')
                curve6.setData(xdata6,pen='b')
                curve6_1.setData(xdata6_1,pen='m')
                curve7.setData(xdata7,pen='r')
                curve8.setData(xdata8,pen='w')
                curve8_1.setData(xdata8_1,pen='m')
                curve9.setData(xdata9,pen='b')
                curve9_1.setData(xdata9_1,pen='m')
                curve10.setData(xdata10,pen='r')
                curve11.setData(xdata11,pen='r')
                ok2=0;

try:
    timer1 = QtCore.QTimer()
    timer1.timeout.connect(update)
    timer1.start(0)
    
    timer2 = QtCore.QTimer()
    timer2.timeout.connect(update2)
    timer2.start(0)

    timer3 = QtCore.QTimer()
    timer3.timeout.connect(Serial)
    timer3.start(0)

    timer4 = QtCore.QTimer()
    timer4.timeout.connect(plots)
    timer4.start(20)
    
except KeyboardInterrupt:
      print('exiting')
      timer.stop()
      arquivo.close();
      raw.close();
    

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
        raw.write("p");
        raw.close();
        arquivo.close();
        arquivo2.close();
        
