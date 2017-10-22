#!/usr/bin/python

import sys
import os
import subprocess
from subprocess import call
import fileinput
import math
import time

global destDir
global imageCounter
global plotCounter

imageCounter = 0
plotCounter = 0
destDir = "compiled/"

def nsplit(string, char, stripp):   #function to split strings preserving quoted passages
    a = []
    out = []
    goOn = True
    if '"' in string:           #split single or double quotes
        a = string.split('"')
    elif "'" in string:
        a = string.split("'")
    else:
        out = string.split(char)  #if no quotes normal split
        return out
        goOn = False
        
    if (len(a) % 2) == 0:
        print "Quotation error!\n"
        goOn = False
        
    if goOn:
        for i in range(0, len(a)-1)[::2]:   #use every other
            b = a[i].split(char)
            if i != 0:
                e = a[i-1] + b[0]       #merge previous quoted passage with first element of unquoted and add to list
                if stripp:
                    out.append(e.strip())
                else:
                    out.append(e)
            if(len(b) > 2):
                out.extend(b[1:len(b)-2])   #add unquoted split passages to list
               
            e = b[len(b)-1] + a[i+1]    #modify next element by adding previous unquoted
            if stripp:
                a[i+1] = e.strip()
            else:
                a[i+1] = e
            
        out.append(a[len(a) - 2])   #add last element
        return out

def makeTable(destFile, args):
    delim = ""
    setDelim = False
    try:
        tabFile = open(args[0])     #add option to pass delimiter in args
        createTableHeader = True
        createTableFooter = False
        
        mergeSet = {}
        mergeNumber = 0
        mathSet = []
        mathAll = False
        caption = ""
        label = ""
        
        for a in args:
            if "merge" in a:
                mArgs = a.split("[")[1].strip("]").split(",")
                try:
                    mergeSet[int(mArgs[0].strip())] = [int(mArgs[1].strip()), mArgs[2]]
                    mergeNumber += (int(mArgs[1].strip()) - int(mArgs[0].strip()))
                except:
                    print "Wrong arguments for 'merge()' in table "+destFile
            if "del" in a:
                delim = a.split("[")[1].strip("]").strip().strip("'").strip('"')
                setDelim = True
            if "math" in a:
                maArgs = a.split("[")[1].strip("]").replace(" ", "").split(",")
                if "all" in maArgs:
                    mathAll = True
                else:
                    try:
                        for arg in maArgs:
                            if "-" in arg:
                                mRange = arg.split("-")
                                for t in range(int(mRange[0]), int(mRange[1])+1):
                                    mathSet.append(t)
                            else:
                                mathSet.append(int(arg.strip()))
                    except ValueError:
                        print "Wrong arguments for 'math()' in table "+destFile
            if "caption" in a:
                caption = a.split("[")[1].strip("]").strip().strip("'").strip('"')
            if "label" in a:
                label = a.split("[")[1].strip("]").strip().strip("'").strip('"')
            if "noheader" in a:
                createTableHeader = False
                    
        for line in tabFile:
            if not setDelim:
                if "\t" in line:
                    delim = "\t"
                elif " " in line:
                    delim = " "
                elif ";" in line:
                    delim = ";"
                elif "|" in line:
                    delim = "|"
                elif "," in line:
                    delim = ","
                elif "." in line:
                    delim = "."
                
                setDelim = True
                print "Delimiter is: "+delim
                
            tabRow = []
            if ('"' in line) or ("'" in line):  #merge strings together
                tempTabRow = line.split(delim)
                i = -1
                a = b = -1
                q = ""
                merge = False
                for e in tempTabRow:
                    i += 1
                    if not merge:
                        if (e.strip()[0] == '"') or (e.strip()[0] == "'"):
                            q = e.strip()[0]
                            if e.strip()[len(e.strip())-1] == q:
                                tabRow.append(e.strip())
                            else:
                                a = i
                                merge = True
                        else:
                            tabRow.append(e.strip())
                    else:
                        if e.strip()[len(e.strip())-1] == q:
                            b = i+1
                            merge = False
                            tabRow.append(delim.join(tempTabRow[a:b]).strip().strip(q))
            else:
                tempTabRow = line.split(delim)
                for e in tempTabRow:
                    tabRow.append(e.strip())
            
            if createTableHeader:           #amount of collumns and delimiter based on first line!  Add option noheader, so only table cells are converted
                createTableHeader = False
                createTableFooter = True
                destFile.write("\\begin{table}[!hbt]\n")
                destFile.write("\\center\n\\caption{"+caption+"}\n\\label{"+label+"}\n")  #caption and label
                cols = "|"
                for j in range(len(tabRow)-mergeNumber):
                    cols += "l|"       #add param for custom align - for now just left 
                destFile.write("\\begin{tabular}{"+cols+"}\n")
                destFile.write("\\hline\n")
                destRow = ""
                c = 0
                mergeCols = 0
                for e in tabRow:
                    c+=1
                    if mergeCols > 0:       #Take care of merging cols
                        if mergeSet[mergeCols][0] == c:
                            if e.strip().strip("\t").strip() != "":
                                destRow += e + " & "
                            else:
                                destRow += " & "
                            mergeCols = 0
                        else:
                            if c == len(tabRow):    #if row shortened
                                destRow += e
                            else:
                                if tabRow[c] == "":     #if last cell empty
                                    destRow += e
                                else:
                                    destRow += e + mergeSet[mergeCols][1].strip().strip("'").strip('"')
                    else:
                        if c in mergeSet:
                            if c < mergeSet[c][0]:
                                mergeCols = c
                                if c == len(tabRow):
                                    destRow += e
                                else:
                                    if tabRow[c] == "":
                                        destRow += e
                                    else:
                                        destRow += e + mergeSet[mergeCols][1].strip().strip("'").strip('"')
                            else:
                                destRow += e + " & "
                        else:
                            destRow += e + " & "
                destFile.write(destRow.strip(" & "))    
                destFile.write("\\\\\n")
                destFile.write("\\hline\n")
            else:
                destRow = ""
                c = 0
                mergeCols = 0
                mathCell = False
                for e in tabRow:
                    c+=1
                    if mergeCols > 0:       #Take care of merging cols
                        if mergeSet[mergeCols][0] == c:
                            if mathCell:
                                destRow += e + "$ & "
                                mathCell = False
                            else:
                                destRow += e + " & "
                            mergeCols = 0
                        else:
                            if c == len(tabRow):
                                destRow += e
                            else:
                                if tabRow[c] == "":
                                    destRow += e
                                else:
                                    destRow += e + mergeSet[mergeCols][1].strip().strip("'").strip('"')
                    else:
                        if c in mergeSet:
                            if c < mergeSet[c][0]:
                                mergeCols = c
                                if mathAll or (c in mathSet):
                                    if c == len(tabRow):
                                       destRow += "$"+ e + "$"
                                    else:
                                       if tabRow[c] == "":
                                           destRow += "$" + e
                                       else:
                                           destRow += "$"+ e + mergeSet[mergeCols][1].strip().strip("'").strip('"')
                                    mathCell = True
                                else:
                                    if c == len(tabRow):
                                        destRow += e
                                    else:
                                        if tabRow[c] == "":
                                            destRow += e
                                        else:
                                            destRow += e + mergeSet[mergeCols][1].strip().strip("'").strip('"')
                            else:
                                if mathAll or (c in mathSet):
                                    destRow += "$"+ e + "$ & "
                                else:
                                    destRow += e + " & "
                        else:
                            if mathAll or (c in mathSet):
                                destRow += "$"+ e + "$ & "
                            else:
                                destRow += e + " & "
                                
                destFile.write(destRow.strip(" & "))
                destFile.write("\\\\\n")
        if createTableFooter:
            destFile.write("\\hline\n\\end{tabular}\n\\end{table}\n")
            
        tabFile.close()
    except:
        print "Error opening table "+args[0]
                

def image(destfile, args):
    global imageCounter
    filenames = args[0].strip("[").strip("]").split(",")
    
    caption = []
    label = []
    size = 0.8
    maxOnPage = -1
    nobreak = False
    doEps = False
    autolabel = False
    
    for a in args:
        if "caption" in a:
            caption = nsplit(a.split("[")[1].strip("]"), ",", True)
        if "label" in a and not "autolabel" in a:
            label = nsplit(a.split("[")[1].strip("]"), ",", True)
        if "size" in a:
            try:
                size = float(a.split("[")[1].strip("]").strip().strip("'").strip('"').replace(" ", ""))    #must be number to avoid latex errors
            except ValueError:
                size = 0.8        #default size
        if "nobreak" in a:
            nobreak = True
        if "autolabel" in a:
            autolabel = True
        if "maxonpage" in a:
            try:
                maxOnPage = int(a.split("[")[1].strip("]").strip().strip("'").strip('"').replace(" ", ""))    #must be number to avoid errors
            except ValueError:
                pass
    
    try:
        if not os.path.exists(destDir+"images"):
            call(["mkdir", destDir+"images"])
    except:
        print "Error creating image directory"
    if len(filenames) == 1:
        filename = filenames[0].strip().strip("'").strip('"')
        try:
            if os.path.isfile(filename):
                imageCounter += 1
                doEps = False
                fileExt = filename.split(".")[len(filename.split("."))-1]
                if fileExt != "png" and fileExt != "eps":
                    call(["convert", filename, destDir+"images/"+str(imageCounter)+".png"])       #include resizing, turning etc
                else:
                    call(["/bin/cp", "-f", filename, destDir+"images/"+str(imageCounter)+"."+fileExt])
                    if fileExt == "eps":
                        doEps = True
                
                destfile.write("\\begin{figure}[!hbt]\n\\begin{center}\n")
                if doEps:
                    destfile.write("\\resizebox{"+str(size)+"\\linewidth}{!}{\\includegraphics[1.0\\linewidth]{"+destDir+"images/"+str(imageCounter)+".eps}}\n")
                else:
                    destfile.write("\\includegraphics[width="+str(size)+"\\linewidth]{"+destDir+"images/"+str(imageCounter)+".png}\n")
                if len(caption) > 0:
                    destfile.write("\\caption{"+caption[0].strip().strip("'").strip('"')+"}\n")
                if len(label) > 0:
                    destfile.write("\\label{"+label[0].strip().strip("'").strip('"')+"}\n")
                destfile.write("\\end{center}\n\\end{figure}\n")
            else:
                print "Warning: File "+filename+" not found!"
        except:
            print "Error opening file "+filename
    
    elif len(filenames) > 1:
        if size > 0.5:                  #if size larger 0.5: 2 per page, centering
            centering = "\\centering\n"
            imagesize = "1.0"
            a = 2
        else:
            a = 8
            imagesize = "1.0"
            centering = ""
            
        if maxOnPage != -1 and a > maxOnPage:       #Limit number of images on one page
            a = maxOnPage
            
        c = int(math.ceil(float(len(filenames)) / float(a)))        #a: Anzahl pro Seite, c: Anzahl Seiten

        for i in range(0, c):
            if not nobreak:
                destfile.write("\\clearpage\n")
            destfile.write("\\begin{figure}[!hbt]\n")
            if i == 0:
                if len(label) > 0:
                    destfile.write("%\\label{"+label[0].strip().strip("'").strip('"')+"}\n")
            else:
                destfile.write("\\ContinuedFloat\n")
            
            destfile.write(centering)
            
            for j in range(1, a+1):
                b = i*a + j
                if len(filenames) >= b:
                    filename = filenames[b-1].strip().strip("'").strip('"')
                    try:
                        if os.path.isfile(filename):
                            imageCounter += 1
                            doEps = False
                            fileExt = filename.split(".")[len(filename.split("."))-1]
                            if fileExt != "png" and fileExt != "eps":
                                call(["convert", filename, destDir+"images/"+str(imageCounter)+".png"])       #include resizing, turning etc
                            else:
                                call(["/bin/cp", "-f", filename, destDir+"images/"+str(imageCounter)+"."+fileExt])
                                if fileExt == "eps":
                                    doEps = True
                                    
                            destfile.write("\\begin{subfigure}{"+str(size)+"\\linewidth}\n")
                            
                            if doEps:
                                destfile.write("\\resizebox{"+imagesize+"\\linewidth}{!}{\\includegraphics[1.0\\linewidth]{"+destDir+"images/"+str(imageCounter)+".eps}}\n")
                            else:
                                destfile.write("\\includegraphics[width="+imagesize+"\\linewidth]{"+destDir+"images/"+str(imageCounter)+".png}\n")
                            
                            if len(caption) >= b+1:
                                destfile.write("\\caption{"+caption[b].strip().strip("'").strip('"')+"}\n")
                            elif autolabel and len(label)>0:
                                destfile.write("\\caption{}\n")
                                
                            if len(label) >= b+1:
                                destfile.write("\\label{"+label[b].strip().strip("'").strip('"')+"}\n")
                            elif autolabel and len(label)>0:
                                destfile.write("\\label{"+label[0]+"_"+str(b)+"}\n")
                                
                            destfile.write("\\end{subfigure}\n")
                    except:
                        print "Error opening file " +filename
            if len(caption) > 0:
                destfile.write("\\caption{"+caption[0].strip().strip("'").strip('"')+"}\n")
            if i == 0:
                if len(label) > 0:
                    destfile.write("\\label{"+label[0].strip().strip("'").strip('"')+"}\n")
            destfile.write("\\end{figure}\n")


def plot(destfile, args):       #Plot umschreiben, so dass es subfigure nutzt; ggf auch mehrere gpts
    global plotCounter
    outputCounter = 0
    i = 0
    output = 0
    files = args[0].strip("'").strip('"').strip().split(",")
    caption = [""]
    label = [""]
    colour = ""
    size = 1.0
    
    for a in args:
        if "caption" in a:
            caption = nsplit(a.split("[")[1].strip("]"), ",", True)
        if "colour" in a:
        	colour = " color"
    	if "size" in a:
            try:
                size = float(a.split("[")[1].strip("]").strip().strip("'").strip('"').replace(" ", ""))    #must be number to avoid latex errors
            except ValueError:
                size = 1.0        #default size
        if "label" in a:
            label = nsplit(a.split("[")[1].strip("]"), ",", True)
        if "output" in a:
            outArg = a.split("[")[1].strip("]").strip().strip("'").strip('"')
            if outArg == "all":
                output = 0
            else:
                try:
                    output = int(outArg)
                except ValueError:
                    print "Warning: Wrong argument for output in plot "+files[0]
                    output = -1
            
    try:
        if not os.path.exists(destDir+"plots"):
            call(["mkdir", destDir+"plots"])
        
        plotCounter += 1
        for e in files:
            filename = e.strip().strip("'").strip('"')
            i += 1
            if os.path.isfile(filename):
                gnuscript = open(filename)
                new_gnuscript = open(os.path.dirname(filename)+"/.convLatexPlot"+str(plotCounter)+"_"+str(i)+".gpt", "w")
                for line in gnuscript:
                    if (line[0] != '#') and ("set terminal" in line):
                        new_line = line.split("set terminal")[0] + "set terminal epslatex" + colour + ";" + "".join("".join(line.split("set terminal")[1:]).split(";")[1:])
                    else:
                        new_line = line
                    if (line[0] != '#') and ("set output" in new_line):
                        outputCounter += 1
                        outpath = os.path.realpath("./")+"/"+destDir+"plots/"+str(plotCounter)+"_"+str(outputCounter)+".tex"
                        new_line = new_line.split("set output")[0] + "set output '" + outpath + "';" + "".join("".join(new_line.split("set output")[1:]).split(";")[1:])
                    new_gnuscript.write(new_line)
                new_gnuscript.close()
                gnuscript.close()
                
                scriptDir = os.path.dirname(os.path.realpath(__file__))
                os.chdir(os.path.dirname(filename))
                call(["gnuplot", ".convLatexPlot"+str(plotCounter)+"_"+str(i)+".gpt"])
                call(["rm", ".convLatexPlot"+str(plotCounter)+"_"+str(i)+".gpt"])
                os.chdir(scriptDir)
            else:
                print "Warning: File not found: "+filename
            
        if output == 0:
            if outputCounter > 1:
                destfile.write("\\clearpage\n\\begin{figure}[!hbt]\n\\caption{"+caption[0].strip().strip("'").strip('"')+"}\n\\label{"+label[0].strip().strip("'").strip('"')+"}\n")
                if (size >= 0.5):
                    destfile.write("\\centering\n")
                for i in range(1, outputCounter+1):
                    destfile.write("\\begin{subfigure}{"+str(size)+"\\linewidth}\n\\resizebox{1.0\\linewidth}{!}{\\input{"+destDir+"plots/"+str(plotCounter)+"_"+str(i)+".tex}}\n")
                    if i < len(caption):
                        if caption[i] != "":
                            destfile.write("\\caption{"+caption[i].strip().strip("'").strip('"')+"}\n")
                    if i < len(label):
                        if label[i] != "":
                            destfile.write("\\label{"+label[i].strip().strip("'").strip('"')+"}\n")
                    destfile.write("\\end{subfigure}\n")
                    if (size >= 0.5):
                        if (i % 2 == 0 and i < outputCounter):
                            destfile.write("\\end{figure}\n\\clearpage\n\n\\begin{figure}\n\\ContinuedFloat\n\\centering\n")
                    else:
                        if (i % 8 == 0 and i < outputCounter):
                            destfile.write("\\end{figure}\n\\clearpage\n\n\\begin{figure}\n\\ContinuedFloat\n")
                destfile.write("\\end{figure}\n\n")
            elif outputCounter > 0:
                i = 1
                destfile.write("\\begin{figure}[!hbt]\n\\centering\n\\resizebox{"+str(size)+"\\linewidth}{!}{\\input{"+destDir+"plots/"+str(plotCounter)+"_"+str(i)+".tex}}\n")
                if i <= len(caption):
                    if caption[i-1] != "":
                        destfile.write("\\caption{"+caption[i-1].strip().strip("'").strip('"')+"}\n")
                if i <= len(label):
                    if label[i-1] != "":
                        destfile.write("\\label{"+label[i-1].strip().strip("'").strip('"')+"}\n")
                destfile.write("\\end{figure}\n\n")
                		                
        elif (output > 0) and (output <= outputCounter):
            destfile.write("\\begin{figure}[!hbt]\n\\begin{center}\en\\resizebox{"+str(size)+"\\linewidth}{!}{\\input{"+destDir+"plots/"+str(plotCounter)+"_"+str(output)+".tex}}\n")
            if caption[0] != "":
                destfile.write("\\caption{"+caption[0].strip().strip("'").strip('"')+"}\n")
            if label[0] != "":
                destfile.write("\\label{"+label[0].strip().strip("'").strip('"')+"}\n")
            destfile.write("\\end{center}\n\\end{figure}\n\n")
        elif output > outputCounter:
            print "Warning: argument output out of range for plot "+files[0]
    except:
        print "Warning: Error opening file "+files[0]

files = []
mainFn = sys.argv[1]
debug = False
doFeynman = False
doImg = True

#Feynman Diagrams
feynDiag = ""

for arg in sys.argv:
    if arg == "debug":
        debug = True    #don't convert anything so line numbers match
    if arg == "doFeynman":
        doFeynman = True
    if arg == "noImg":
        doImg = False
        
if mainFn != "":
    main = open(mainFn)
else:
    print "No Filename!"
    exit()
for line in main:
    if line.strip()[0:9] == "\include{" or line.strip()[0:9] == "\input{":
        files.append(line.split("{")[1].split("}")[0]+".tex")
        
    if line.strip()[0:15] == "\\begin{fmffile}":                 #Process Feynman Diagrams
        feynDiag = line.strip("\\begin{fmffile}{").split("}")[0].strip()

main.close()
files.append(mainFn)

if len(files) > 0:
    #try:
        try:
            if not os.path.exists(destDir):
                call(["mkdir", destDir])
        except:
            print "Unable to create destination directory! Exciting..."
            exit()
            
        #call(["/bin/cp", "-f", mainFn, destDir])    
        
        for curFileName in files:
            curFile = open(curFileName)
            destFile = open(destDir + curFileName, "w")
            for line in curFile:
                if line[:2] == "%!" and not debug:
                    args = "(".join(line.split("(")[1:]).strip().strip(")").split(";")
                    if line[2:].strip()[:5] == "table":                         #!!!Add: Don't split text in quotes (also for [ - brackets - ])
                        makeTable(destFile, args)
                    if line[2:].strip()[:5] == "image" and doImg:
                        image(destFile, args)
                    if line[2:].strip()[:4] == "plot":
                        plot(destFile, args)
                    
                else:
                    destFile.write(line)
            
            destFile.close()        
            curFile.close()
        print "Compilation successful!"
        
        #time.sleep(5)
        
        try:

            if debug:
                if not os.path.exists("debug"):
                    call(["mkdir", destDir+"debug"])
                call(["pdflatex", "-shell-escape", "-output-directory", destDir+"debug", mainFn])
            else:
                call(["pdflatex", "-shell-escape", "-output-directory", destDir, mainFn])
                if(doFeynman and feynDiag != ""):
                    call(["mpost", destDir + feynDiag])
                    call(["pdflatex", "-shell-escape", "-interaction", "nonstopmode", "-output-directory", destDir, mainFn])
            
            if not debug:
                call(["/bin/cp", "-f", destDir+mainFn.strip(".tex")+".pdf", "./"])

        except:
            print "Warning: Couldn't execute pdflatex!"
    #except:
    #    print "Error opening tex files!"
    #    exit()
else:
    print "No files to be edited"
                
        
#Implement error check like "\end missing in line bla" or "closing parenthesis missing line blubb" to avoid annoying latex error messages
#convert images using convert - perhaps add args like "turn" if image upside down etc.
