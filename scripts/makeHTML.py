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

def makeHTML(title, htmlDir, plot_types, description, IFRAMEA, doYields=True):

    print "\n\n ******** Now making HTML pages ******** \n"
    menu=""
    #plot_types = ["comb","h_dist1","h_dist2","h_dist3","fit","4mu","4e"]

    fileList = {}
    for x in plot_types:
      newDir = htmlDir+"/"+x
      createDir(newDir)

      fname = x+".html"
      imgfile = open(fname,"w")
      imgfile.write('<html><head><title>'+x+'</title></head><body style="background-color: Linen">\n')

      fileList[x] = os.listdir(newDir)

      PDFfileList = [f[:-4] for f in os.listdir(newDir) if f.endswith(".pdf")]

      # print PDFfileList
      # print "\n In Making html:", x
      # print fileList[x]

      count = 1
      mod   = 1
      for pl in sorted(fileList[x]):
        if not pl.endswith(".png"): continue
        mod = count % 2
        if pl[:-4] in PDFfileList:
          # print pl[:-3]

          if mod==1: imgfile.write('<nobr><a href='+x+'/'+pl[:-4]+'.pdf><img src='+x+'/'+pl+' width=49%></a>')
          if mod==0: imgfile.write('      <a href='+x+'/'+pl[:-4]+'.pdf><img src='+x+'/'+pl+' width=49%></a></nobr>\n')
        else:
          if mod==1: imgfile.write('<nobr><img src='+x+'/'+pl+' width=49%>')
          if mod==0: imgfile.write('      <img src='+x+'/'+pl+' width=49%></nobr>\n')
        count+=1
      if mod==0: imgfile.write("")
      if mod==1: imgfile.write("</nobr>")

      imgfile.write("</body></html>")
      imgfile.close()
      os.system("mv "+fname+" "+htmlDir)

      menu = menu+"<li><a href=\""+x+".html\" target=\"iframe_a\">"+x+"</a></li>"

    if doYields:
      menu += "<li><a href=\"yields.html\" target=\"iframe_a\">yields</a></li>"

    today = datetime.date.today()
    print today

    message = '<h2>Comments</h2>'
    message+='<ul>'
    for d in description:
      message += "<li>"+d+" </li>"
    if doYields:
      message += "<li><a href=\"yields_all.twiki\">twiki</a></li>"
      message += "<li><a href=\"yields_all.tex\">tex</a></li>"
    message+='</ul>'

    tempfile = open("../scripts/indextemplate.html","r")
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

    os.system("cp ../scripts/style.css "+htmlDir)
    if doYields:

      os.system("cp yields* "+htmlDir)
      os.system("rm yields*")

    print "\n\n *** End of  making HTML pages - all done *** \n"


if __name__ == "__main__":

    pathBase = "/uscms_data/d2/andreypz/html/zgamma/test"

    plot_types =["a","b","c"]
    table = [[]]

    makeTable(table,"twiki")

    comments = ["1",
                "Blah"]

    makeHTML("h &rarr; dalitz decay plots", pathBase, plot_types, comments, "a")

    print "\n\t\t finita la comedia \n"
