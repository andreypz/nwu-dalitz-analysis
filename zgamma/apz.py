#! /usr/bin/env python
from ROOT import *
import utils as u

def alphaPiZ2(f1, c1, globalCut, path):
  print "study alpha/piz particle"
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


