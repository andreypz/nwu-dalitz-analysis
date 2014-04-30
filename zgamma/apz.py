#! /usr/bin/env python
from ROOT import *
import utils as u

def treeInOne(name, var, r1, r2, nBins, globalCut, t,c1,path):
  var1 = str(r1)
  var2 = str(r2)
  m1 = '0'
  m2 = '2.9'
  mllCut = TCut('(m12>'+m1+') && (m12<'+m2+')')
  binWidth = (float(var2)-float(var1))/float(nBins)
  rangeCut = TCut('('+var+'>'+var1+')&&('+var+'<'+var2+')')

  name1 = name+m1+'_'+m2
  t.Draw(var+'>>'+name1+'('+str(nBins)+','+var1+','+var2+')',  mllCut+globalCut+rangeCut, 'goff')

  m1 = '2.9'
  m2 = '3.3'
  mllCut = TCut('(m12>'+m1+') && (m12<'+m2+')')
  name2 = name+m1+'_'+m2
  t.Draw(var+'>>'+name2+'('+str(nBins)+','+var1+','+var2+')',  mllCut+globalCut+rangeCut, 'goff')

  m1 = '3.3'
  m2 = '20'
  mllCut = TCut('(m12>'+m1+') && (m12<'+m2+')')
  name3 = name+m1+'_'+m2
  t.Draw(var+'>>'+name3+'('+str(nBins)+','+var1+','+var2+')',  mllCut+globalCut+rangeCut, 'goff')

  h1 = gDirectory.Get(name1)
  h2 = gDirectory.Get(name2)
  h3 = gDirectory.Get(name3)
  h1.Draw("hist")
  h2.Draw("sames hist")
  h3.Draw("sames hist")
  mm = max([h1.GetMaximum(), h2.GetMaximum(), h3.GetMaximum()])
  h1.SetMaximum(1.4*mm)
  h1.SetMinimum(0)

  h2.SetLineColor(kGreen+2)
  h3.SetLineColor(kRed+1)
  #h1.UseCurrentStyle()

  gStyle.SetOptStat(1111)
  h1.SetName("Low Mll")
  h2.SetName("J/Psi")
  h3.SetName("High Mll")
  stats1 = h1.GetListOfFunctions().FindObject("stats");
  stats2 = h2.GetListOfFunctions().FindObject("stats");
  stats3 = h3.GetListOfFunctions().FindObject("stats");
  # stats1.Print()
  stats1.SetX1NDC(0.7)
  stats1.SetX2NDC(0.95)
  stats1.SetY1NDC(0.9)
  stats1.SetY2NDC(0.7)
  stats1.SetTextColor(kBlue+2)

  stats2.SetX1NDC(0.7)
  stats2.SetX2NDC(0.95)
  stats2.SetY1NDC(0.7)
  stats2.SetY2NDC(0.5)
  stats2.SetTextColor(kGreen+2)
  stats3.SetX1NDC(0.7)
  stats3.SetX2NDC(0.95)
  stats3.SetY1NDC(0.5)
  stats3.SetY2NDC(0.3)
  stats3.SetTextColor(kRed+1)

  leg = TLegend(0.20,0.76,0.48,0.92);
  leg.AddEntry(h1,'0.2 < m_{#mu#mu} < 2.9', "l")
  leg.AddEntry(h2,'2.9 < m_{#mu#mu} < 3.3', "l")
  leg.AddEntry(h3,'3.3 < m_{#mu#mu} < 20', "l")
  leg.SetTextSize(0.03)
  leg.SetFillColor(kWhite)
  leg.Draw()

  h1.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

def ZJPG(f1, c1, globalCut, path):
  print "\n\n *** Study of Z/H -> J/Psi Gamma and more ***\n"
  print globalCut

  u.createDir(path)
  c1.cd()
  t = f1.Get('apzTree/apzTree')
  gStyle.SetMarkerSize(0.6)
  gStyle.SetMarkerStyle(20);
  gStyle.SetLineWidth(1);
  gStyle.SetOptStat(111)
  opt = ''

  binDownSize = 1

  mllCut = TCut('')

  name = 'h00-mll-HIG14-003'
  m1 = '0'
  m2 = '20'
  mllCut = TCut('m12<'+m2)

  nBins = 50/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m12>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")


  name = 'h01-mll-25'
  m1 = '0'
  m2 = '30'
  mllCut = TCut('m12<'+m2)
  nBins = 100/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m12>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = 'h02-mll-full'
  m1 = '0'
  m2 = '150'
  mllCut = TCut('m12<'+m2)
  nBins = 100/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m12>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = 'h03-mll-jpsi'
  m1 = '2.9'
  m2 = '3.3'
  mllCut = TCut('(m12>'+m1+') && (m12<'+m2+')')
  nBins = 20/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m12>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = "h04-mllg-full"
  mllg1 = '0'
  mllg2 = '200'
  mllCut = TCut('')
  nBins = 50/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  rangeCut = TCut('(m123>'+mllg1+')&&(m123<'+mllg2+')')
  t.Draw('m123>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut+rangeCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h05-mllg-in_Mll30"
  mllg1 = '60'
  mllg2 = '150'
  m1 = '0'
  m2 = '30'
  mllCut = TCut('(m12>'+m1+') && (m12<'+m2+')')
  nBins = 50/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  rangeCut = TCut('(m123>'+mllg1+')&&(m123<'+mllg2+')')
  t.Draw('m123>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut+rangeCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h06-dr12"
  nBins = 50/binDownSize
  dr1 = '0'
  dr2 = '4'
  binWidth = (float(dr2)-float(dr1))/float(nBins)
  rangeCut = TCut('(dr12>'+dr1+')&&(dr12<'+dr2+')')
  t.Draw('dr12>>'+name+'('+str(nBins)+','+dr1+','+dr2+')',  globalCut+rangeCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';#Delta R(#mu_{1}, #mu_{2});Events/%.2f' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h07-dr1234"
  nBins = 50/binDownSize
  dr1 = '0'
  dr2 = '4'
  binWidth = (float(dr2)-float(dr1))/float(nBins)
  rangeCut = TCut('(dr1234>'+dr1+')&&(dr1234<'+dr2+')')
  t.Draw('dr1234>>'+name+'('+str(nBins)+','+dr1+','+dr2+')',  globalCut+rangeCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';#Delta R(#mu#mu, #gamma);Events/%.2f' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  treeInOne('h08-mllg_full_threeInOne', 'm123',  0, 200, 50, globalCut, t,c1,path)
  treeInOne('h09-mllg_mZ_threeInOne',   'm123', 50, 140, 50, globalCut, t,c1,path)
  treeInOne('h10-mllg_Hfull_threeInOne','m123',110, 170, 20, globalCut, t,c1,path)
  treeInOne('h11-mllg_mH_threeInOne',   'm123',115, 135, 10, globalCut, t,c1,path)

def alphaPiZ2(f1, c1, globalCut, path):
  print "\n\n *** Study alpha/piz particle ***\n"
  u.createDir(path)
  c1.cd()
  t = f1.Get('apzTree/apzTree')
  gStyle.SetMarkerSize(0.6)
  gStyle.SetMarkerStyle(20);
  gStyle.SetLineWidth(1);
  gStyle.SetOptStat(111)
  opt = ''

  binDownSize = 1

  name = 'h01-m12-full'
  mllCut = TCut('')
  # myFullCut = TCut(mllCut+globalCut)
  m1 = '0'
  m2 = '150'
  nBins = 100/binDownSize
  # print 10*'***', 'Downsized to ', nBins
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m12>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{12} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")


  name = 'h02-m34-full'
  mllCut = TCut('')
  # myFullCut = TCut(mllCut+globalCut)
  m1 = '0'
  m2 = '150'
  nBins = 100/binDownSize
  # print 10*'***', 'Downsized to ', nBins
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m34>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{34} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")


  name = 'h03-ml12-apz'
  mllCut = TCut('')
  # myFullCut = TCut(mllCut+globalCut)
  m1 = '10'
  m2 = '35'
  nBins = 15/binDownSize
  # print 10*'***', 'Downsized to ', nBins
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m12>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{12} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = 'h04-ml34-apz'
  mllCut = TCut('')
  # myFullCut = TCut(mllCut+globalCut)
  m1 = '10'
  m2 = '35'
  nBins = 15/binDownSize
  # print 10*'***', 'Downsized to ', nBins
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m34>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{12} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")


  name = 'h05-m4l_full'
  mllCut = TCut('')
  m1 = '60'
  m2 = '200'
  nBins = 30/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m4l>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{4l} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")


  name = 'h06-m4l_full'
  mllCut = TCut('')
  m1 = '100'
  m2 = '180'
  nBins = 20/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m4l>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{4l} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")



  name = 'h07-pt12'
  mllCut = TCut('')
  m1 = '0'
  m2 = '180'
  nBins = 50/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('pt12>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';p_{T}^{12} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = 'h08-pt34'
  mllCut = TCut('')
  m1 = '0'
  m2 = '180'
  nBins = 50/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('pt34>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';p_{T}^{34} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = 'h09-pt12_over_m4l'
  mllCut = TCut('')
  m1 = '0'
  m2 = '1'
  nBins = 50/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('pt12/m4l>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';p_{T}^{12}/m_{4l} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = 'h10-pt34_over_m4l'
  mllCut = TCut('')
  m1 = '0'
  m2 = '1'
  nBins = 50/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('pt34/m4l>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';p_{T}^{34}/m_{4l} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")



def alphaPiZ(f1, c1, globalCut, path):
  print "study alpha/piz particle"
  u.createDir(path)
  c1.cd()
  t = f1.Get('fitTree/fitTree')
  gStyle.SetMarkerSize(0.6)
  gStyle.SetMarkerStyle(20);
  gStyle.SetLineWidth(1);
  gStyle.SetOptStat(111)
  opt = ''

  binDownSize = 1
  if 'alphaPiZ-3' in path:
    binDownSize = 2
  elif 'alphaPiZ-4' in path:
    binDownSize = 2
  elif 'alphaPiZ-6' in path:
    binDownSize = 2

  name = 'h00-mll-HIG14-003'
  mllCut = TCut('')
  # myFullCut = TCut(mllCut+globalCut)
  m1 = '0'
  m2 = '20'
  nBins = 50/binDownSize
  # print 10*'***', 'Downsized to ', nBins
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m_ll>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")



  name = 'h01-mll-full'
  mllCut = TCut('')
  # myFullCut = TCut(mllCut+globalCut)
  m1 = '0'
  m2 = '25'
  nBins = 100/binDownSize
  # print 10*'***', 'Downsized to ', nBins
  binWidth = (float(m2)-float(m1))/float(nBins)
  t.Draw('m_ll>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  h.UseCurrentStyle()
  c1.SaveAs(path+"/"+name+".png")

  name = 'h02-mll-jpsi'
  m1 = '2.7'
  m2 = '4.0'
  nBins = 100/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_ll>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h03-mll-upsilon"
  m1 = '9.0'
  m2 = '12.0'
  nBins = 30/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_ll>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h04-mll-apz"
  m1 = '15.0'
  m2 = '25.0'
  nBins = 30/binDownSize
  binWidth = (float(m2)-float(m1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_ll>>'+name+'('+str(nBins)+','+m1+','+m2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h05-mllg-full"
  mllg1 = '60'
  mllg2 = '150'
  m1 = '0'
  m2 = '25'
  nBins = 50/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h06-mllg-jpsi-full"
  mllg1 = '100'
  mllg2 = '150'
  m1 = '3.0'
  m2 = '3.2'
  nBins = 20/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('m_ll>'+m1+' && m_ll<'+m2)
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',mllCut*globalCut, opt)
  # print 3*"***", 'my integrall =', h.Integral()
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

  name = "h07-mllg-jpsi-Z"
  mllg1 = '70'
  mllg2 = '110'
  m1 = '3.0'
  m2 = '3.2'
  nBins = 20/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h08-mllg-upsilon-Z"
  mllg1 = '70'
  mllg2 = '110'
  m1 = '9.2'
  m2 = '9.6'
  nBins = '30'
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

  name = "h09-mllg-apz-full"
  mllg1 = '70'
  mllg2 = '200'
  m1 = '17.0'
  m2 = '19.0'
  nBins = 100/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")



  name = "h10-mllg-apz-full"
  mllg1 = '100'
  mllg2 = '150'
  m1 = '17.0'
  m2 = '19.0'
  nBins = 20/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('m_ll>'+m1+' && m_ll<'+m2)
  # print 15*"*", 'myCut is = ', myFullCut
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  # t.Draw('m_llg>>'+name,  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h11-mllg-apz-Z"
  mllg1 = '70'
  mllg2 = '110'
  m1 = '17'
  m2 = '19'
  nBins = 40/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  myFullCut = TCut(mllCut+globalCut)
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


  name = "h12-phEta-apz"
  m1 = '17.0'
  m2 = '19.0'
  nBins = 50/binDownSize
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  myFullCut = TCut(mllCut+globalCut)
  t.Draw('fabs(ph_eta)>>'+name+'('+str(nBins)+',0,3)',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';|#eta_{#gamma}|;Events')
  c1.SaveAs(path+"/"+name+".png")

  name = "h13-mllg-apz-125"
  mllg1 = '115'
  mllg2 = '135'
  m1 = '17.0'
  m2 = '19.0'
  nBins = 30/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  myFullCut = TCut(mllCut+globalCut)
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

  name = "h14-mllg-jpsi-125"
  mllg1 = '115'
  mllg2 = '135'
  myFullCut = TCut(mllCut+globalCut)
  m1 = '3.0'
  m2 = '3.2'
  nBins = 30/binDownSize
  binWidth = (float(mllg2)-float(mllg1))/float(nBins)
  mllCut = TCut('(m_ll>'+m1+') && (m_ll<'+m2+')')
  t.Draw('m_llg>>'+name+'('+str(nBins)+','+mllg1+','+mllg2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';m_{#mu#mu#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

  name = "h15-ptgamma-apz"
  pt1='40'
  pt2='120'
  m1 = '17.0'
  m2 = '19.0'
  nBins = 30/binDownSize
  binWidth = (float(pt2)-float(pt1))/float(nBins)
  mllCut = TCut('m_ll>'+m1+' && m_ll<'+m2)
  myFullCut = TCut(mllCut+globalCut)
  t.Draw('ph_pt>>'+name+'('+str(nBins)+','+pt1+','+pt2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';p_{T}^{#gamma} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

  name = "h16-ptll-apz"
  pt1='40'
  pt2='120'
  m1 = '17.0'
  m2 = '19.0'
  nBins = 30/binDownSize
  binWidth = (float(pt2)-float(pt1))/float(nBins)
  mllCut = TCut('m_ll>'+m1+' && m_ll<'+m2)
  myFullCut = TCut(mllCut+globalCut)
  t.Draw('di_pt>>'+name+'('+str(nBins)+','+pt1+','+pt2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';p_{T}^{#mu#mu} (GeV);Events/%.2f GeV' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

  name = "h17-ptgamma-overMllg-apz"
  pt1='0.3'
  pt2='1.0'
  m1 = '17.0'
  m2 = '19.0'
  nBins = 30/binDownSize
  binWidth = (float(pt2)-float(pt1))/float(nBins)
  mllCut = TCut('m_ll>'+m1+' && m_ll<'+m2)
  myFullCut = TCut(mllCut+globalCut)
  t.Draw('ph_pt/m_llg>>'+name+'('+str(nBins)+','+pt1+','+pt2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';p_{T}^{#gamma}/m_{#mu#mu#gamma};Events/%.2f bin' % binWidth)
  c1.SaveAs(path+"/"+name+".png")

  name = "h18-ptll-overMllg-apz"
  pt1='0.3'
  pt2='1.0'
  m1 = '17.0'
  m2 = '19.0'
  nBins = 30/binDownSize
  binWidth = (float(pt2)-float(pt1))/float(nBins)
  mllCut = TCut('m_ll>'+m1+' && m_ll<'+m2)
  myFullCut = TCut(mllCut+globalCut)
  t.Draw('di_pt/m_llg>>'+name+'('+str(nBins)+','+pt1+','+pt2+')',  mllCut+globalCut, opt)
  h = gDirectory.Get(name)
  h.Draw("same e1p")
  h.UseCurrentStyle()
  h.SetTitle(name+';p_{T}^{#mu#mu}/m_{#mu#mu#gamma};Events/%.2f bin' % binWidth)
  c1.SaveAs(path+"/"+name+".png")


