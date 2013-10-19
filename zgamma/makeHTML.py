#! /usr/bin/env python
import sys,os,datetime

def createDir(dir):
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except OSError:
            if os.path.isdir(dir):
                pass
            else:
                raise
                                            
def makeHTML(title, htmlDir, plot_types, description, IFRAMEA):

    print "\n\n ******** Now making HTML pages ******** \n"
    menu=""
    #plot_types = ["comb","h_dist1","h_dist2","h_dist3","fit","4mu","4e"]

    fileList = {}
    for x in plot_types:
        newDir = htmlDir+"/"+x
        createDir(newDir)
        #if not os.path.exists(newDir):
        #    print "Creating ", newDir
        #    os.makedirs(newDir)


        fname = x+".html"
        imgfile = open(fname,"w")
        imgfile.write("<html><head><title>"+x+"</title></head><body>\n")

        fileList[x] = os.listdir(newDir)
        #print fileList[x]
        
        count =1
        mod =1
        for pl in  sorted(fileList[x]):
            mod = count % 2
            if mod==1: imgfile.write('<nobr><img src='+x+'/'+pl+' width=45%>')
            if mod==0: imgfile.write('      <img src='+x+'/'+pl+' width=45%></nobr>\n')
            count+=1
        if mod==0: imgfile.write("")
        if mod==1: imgfile.write("</nobr>")

        imgfile.write("</body></html>")
        imgfile.close()
        os.system("mv "+fname+" "+htmlDir)

        menu = menu+"<li><a href=\""+x+".html\" target=\"iframe_a\">"+x+"</a></li>"

    menu += "<li><a href=\"yields.html\" target=\"iframe_a\">yields</a></li>"

    today = datetime.date.today()
    print today

    message = '<h2>Comments</h2>'
    message+='<ul>'
    for d in description:
        message+="<li>"+d+" </li>"
        message+='</ul>'

    tempfile = open("indextemplate.html","r")
    whole_thing = tempfile.read()
    whole_thing = whole_thing.replace("{TITLE}", title)
    whole_thing = whole_thing.replace("{MENU}", menu)
    whole_thing = whole_thing.replace("{DATE}", str(today))
    whole_thing = whole_thing.replace("{ASIDEMESSAGE}", str(message))
    whole_thing = whole_thing.replace("{IFRAMEA}", IFRAMEA)
    tempfile.close()


    ifile = open(htmlDir+"/index.html","w")
    ifile.write(whole_thing)
    ifile.close()

    os.system("cp yields.html "+htmlDir)
    
    print "\n\n *** End of  making HTML pages - all done *** \n"


def makeTable(table, opt="tex"):
    print "Making sure that the list is alright"
    n_row = len(table)
    n_col = len(table[0])
    for l in table:
        if len(l)!=n_col:
            print "No good, the number of columns is messed up"
            
            
    myTable = ''
    if opt=="tex":
        beginTable = '\\begin{tabular}{|'+n_col*"l|"+'} \n \\hline \n'
        endTable   = '\\end{tabular} \n'

        beginLine  = ''
        endLine    = '\\\\\\hline \n'
        separator  = ' & '


    if opt=="html":
        beginTable = '<table border = "10"    cellpadding="5">'
        endTable = '</table>'

        beginLine  = '\n<tr>\n<td>'
        endLine    = '</td>\n</tr>'
        separator  = '</td><td>'

    if opt=="twiki":
        beginTable = ''
        endTable   = ''

        beginLine  = '| '
        endLine    = ' |\n'
        separator  = ' |  '


    myTable +=beginTable
    for l in range(n_row):
        #print l, table[l]
        myTable+=beginLine
        for c in range(n_col):
            val = table[l][c]
            if not isinstance(val,str):
                myTable+="%.0f" % (table[l][c])
            else:
                myTable+=val
            if c!=n_col-1:
                myTable+=separator
        myTable+=endLine

    myTable +=endTable
    
    ifile = open("yields.html","w")
    ifile.write(myTable)
    ifile.close()
    print myTable


if __name__ == "__main__":
    
    pathBase = "/uscms_data/d2/andreypz/html/zgamma/test"

    plot_types =["a","b","c"]
    table = [[]]
        
    makeTable(table,"twiki")

    comments = ["1",
                "Blah"]
    
    makeHTML("h &rarr; dalitz decay plots",pathBase, plot_types, comments, "a")

    print "\n\t\t finita la comedia \n"
