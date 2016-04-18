#!/usr/bin/env python
# usage: python convertFits2ROOT.py <inputfile> 

import ROOT as Rt
import numpy as np
import pyfits
import sys
import glob
import array 
import copy, os
if not sys.flags.interactive:
    #print 'batch mode'
    Rt.gROOT.SetBatch(True)
if not Rt.gROOT.GetStyle("Modern"): # checks if Modern is installed (past 5.30)
    Rt.gROOT.SetStyle("Plain") # fall back
else:
    Rt.gROOT.SetStyle("Modern") # otherwise use it!

bold = 1 # this is used to make lines thicker...

def ContainsExpr(string,expr):
    ''' returns True if the expression expr is contained in string '''
    if isinstance(expr,list):
        istrue = False
        for e in expr:
            if e in string:
                istrue = True
        return istrue
    try:
        parts = string.partition(expr);
        if (len(parts[1])>0):
            return True;
        else:
            return False;
    except IndexError:
        False;

def main(fitsfile,rootfile = None,img_extension='png',xpix=600,ypix=600,iTitle=None,catalog_names="2FGL",show_counts=False,rebin=None,plot_total=False):
    cat_names = catalog_names.split(";")
    Rt.gStyle.SetOptStat(0)
    inputfile = fitsfile
    hdu = pyfits.open(inputfile)
    

    spectra = hdu[1]
    corner_title = None
    if dict(hdu[0].header).has_key("ROINAME"):
        corner_title = int(hdu[0].header.get("ROINAME"))

    ebounds = hdu[3].data[:]
    E_MIN = list(hdu[3].data.field("E_MIN"))
    E_MAX = list(hdu[3].data.field("E_MAX"))

    E_MIN.append(E_MAX[-1])
    #print E_MIN 
    #print E_MAX


    hists = []
    nHist = len(hdu[1].data[0])
    for i in range(nHist):
         title = hdu[1].header['TTYPE%s'%str(i+1)]
         name = title.replace(" ","_").replace("-","m").replace(".","d").replace("+","p").replace("2FGL","SRC_2FGL")
         hists.append(Rt.TH1D(name,title,len(E_MIN)-1,array.array("d",E_MIN)))

    #h = hists[0]
    # #print hdu[1].data

    for h in hists:
        dataset = hdu[1].data.field(h.GetTitle())
        #print dataset
        for i in range(len(dataset)):
            h.SetBinContent(i+1,dataset[i])

    htot = Rt.TH1D("Total","Sum of Models",len(E_MIN)-1,array.array("d",E_MIN))
    htot.SetLineColor(1)
    diff = Rt.TH1D("residual","(counts-model)/model",len(E_MIN)-1,array.array("d",E_MIN))
    f1 = Rt.TF1("zero","[0]",E_MIN[0],E_MIN[-1])
    f1.SetParameter(0,0)
    f1.SetLineStyle(3)
    for i in range(len(hists)-1):
        #print hists[i+1].GetTitle()
        htot.Add(hists[i+1])

    htot.SetLineStyle(2)
    if show_counts:
        print '*INFO* Integral over model counts %s : %1.2f'%(fitsfile,float(htot.Integral()))

    Min = hists[0].GetMinimum()
    Max = hists[0].GetMaximum()
    if Min == 0:
        Min = 1e-2
    htot.GetYaxis().SetRangeUser(0.8*Min,1.5*Max)#8e-3,8)
    htot.GetXaxis().CenterTitle()
    htot.GetYaxis().CenterTitle()
    htot.GetXaxis().SetTitle("Energy (MeV)")
    htot.GetYaxis().SetTitle("counts/bin")
    htot.GetYaxis().SetTitleSize(0.06)
    htot.GetYaxis().SetTitleOffset(0.89)
    htot.GetYaxis().SetLabelSize(0.05)
    diff.GetXaxis().SetTitle("Energy (MeV)")
    diff.GetYaxis().SetTitle("#frac{counts-model}{model}")
    diff.GetXaxis().CenterTitle()
    diff.GetYaxis().CenterTitle()

    diff.GetXaxis().SetTitleSize(0.1)
    diff.GetYaxis().SetTitleSize(0.1)
    diff.GetXaxis().SetTitleSize(0.1)
    diff.GetYaxis().SetTitleOffset(0.5)
    diff.GetXaxis().SetTitleOffset(0.95)
    diff.GetYaxis().SetLabelOffset(0.012)
    diff.GetXaxis().SetTickLength(0.09)
    diff.GetYaxis().SetLabelSize(0.10)
    diff.GetXaxis().SetLabelSize(0.08)
    diff.SetLineColor(1)

    counts = hists[0]
    counts.SetMarkerStyle(2)
    #counts.SetMarkerSize(10)
    #counts.SetMarkerColor(1)

    diff.Add(counts)
    diff.Add(htot,-1)
    htot2 = copy.copy(htot)
    htot2.Multiply(htot)
    diff.Divide(htot)
    htot2.Sumw2()
    diff.Sumw2()
    #Rt.gStyle.SetPadTickX(1)
    #Rt.gStyle.SetPadTickY(1)

    #residual = counts_model.Divide(htot)
    #http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=6980
    c1 = Rt.TCanvas("c1","c1",int(xpix),int(ypix))
    p1 = Rt.TPad("pad1","pad1",0,0.33,1,1)
    p2 = Rt.TPad("pad2","pad2",0,0,1,0.33)
    p1.SetLineWidth(int(p1.GetLineWidth()+bold))
    p2.SetLineWidth(int(p2.GetLineWidth()+bold))
    p1.SetBottomMargin(0.00001)
    p1.SetLeftMargin(0.11)
    p1.SetBorderMode(0)
    p1.SetLogy()
    p1.SetLogx()

    p2.SetTopMargin(0)
    p2.SetBottomMargin(0.2)
    p2.SetLeftMargin(0.11)
    p2.SetBorderMode(0)
    p1.Draw()
    p2.Draw()

    #c1 = Rt.TCanvas("c1","canvas",400,600)
    #c1.Divide(1,2)#,0.01,0)
    p1.cd()
    p1.SetLogx()
    p1.SetLogy()
    Leg = Rt.TLegend(0.53,0.885,0.875,0.69)
    Leg.SetFillColor(0)
    Leg.SetLineColor(0)
    Leg.SetTextSize(0.04)
    #htot.GetYaxis().SetRangeUser(
    Leg.AddEntry(counts,'Counts','P')
    Leg.AddEntry(htot,'Sum of Models','l')
    htot.SetLineWidth(int(htot.GetLineWidth()+bold))
    if iTitle is None:
        htot.SetTitle("")
    elif iTitle=="title":
        htot.SetTitle(os.path.basename(fitsfile))
    else:
        #print '*DEBUG* use this title %s'%iTitle
        htot.SetTitle(str(iTitle))
    # if not corner_title is None:
    #     #htot.SetTitleOffset(-0.2)
    #     htot.SetTitle("%i"%corner_title)

    p1.SetTicky()
    p2.SetTicky()
    #p2.SetTickx()
    #htot.GetYaxis().SetNdivisions(710)
    htot.Draw("L")
    added = False
    added_GCl = False
    for i in range(len(hists)):
        hists[i].SetLineWidth(int(hists[i].GetLineWidth()+bold))
        drawoption = "Lsame"
        if i==0:
            drawoption = "PE1same"
        if hists[i].GetTitle()=="EGAL":
             Leg.AddEntry(hists[i],'Extragalactic Diffuse','l')
             hists[i].SetLineColor(2)
             hists[i].SetLineStyle(2)
        elif hists[i].GetTitle()=="GAL":
             Leg.AddEntry(hists[i],'Galactic Diffuse','l')
             hists[i].SetLineColor(3)
             hists[i].SetLineStyle(4)
        elif not ContainsExpr(hists[i].GetTitle(),cat_names) and not hists[i].GetTitle()=="ObsCounts":
            
             hists[i].SetLineColor(4)
             hists[i].SetLineStyle(4)
             if not added_GCl:
                  #Leg.AddEntry(hists[i],'Extra Component','l')
                  added_GCl = True
        else:
             if not added:
                  Leg.AddEntry(hists[i],'Individual Source Models','l')
                  added = True
             hists[i].SetLineColor(1)
        if plot_total and i>0:
            continue
        hists[i].Draw(drawoption)
    #htot.Sumw2()
    hists.append(htot)
    #Leg.Draw("same")
    if not corner_title is None:
        tl = Rt.TLatex()
        tl.SetTextFont(62)
        tl.SetNDC()
        #tl.DrawLatex(0.5,0.91,"%i"%corner_title)
        tl.DrawLatex(0.8,0.8,"%i"%corner_title)
    p2.cd()
    p2.SetLogx()
    #p2.SetCanvasSize(p2.GetWw(),int(p2.GetWh()*0.25))
    # diff.GetYaxis().SetLimits(-0.75,0.75)
    # diff.GetYaxis().SetRangeUser(-0.75,0.75)
    diff.GetYaxis().SetLimits(-1.25,1.25)
    diff.GetYaxis().SetRangeUser(-1.25,1.25)
    diff.GetYaxis().SetNdivisions(206)
    diff.SetLineWidth(int(diff.GetLineWidth()+bold))
    diff.SetTitle("")
    diff.Draw("PE1")
    f1.SetLineWidth(int(f1.GetLineWidth()+bold))
    f1.Draw("lsame")

    if not rootfile is None:
        rh = Rt.TFile(rootfile,'recreate')
        for h in hists:
            h.Write()
        diff.Write()
        c1.Write()
        rh.Write()
        rh.Close()
        print '*wrote ROOT file: %s'%rootfile
    c1.SaveAs(inputfile.replace(".fits","."+img_extension))
    return c1
    
if __name__ == "__main__":
    """ code to convert ScienceTools counts_spectra.fits files to ROOT and pretty plots """
    from optparse import OptionParser    
    parser = OptionParser()
    parser.add_option("--rootfile",dest="rootfile",type=str, default = None, help="name of the rootfile to be written, if run on mutliple files,create a rootfile with the same name as the input but root extension")
    parser.add_option("--img_extension",dest="img_extension",type=str,default = "png", help="image extension, supports png, eps, ps, pdf, jpg?")
    parser.add_option("--x_pix",dest="xpix", type=int, default = 600, help="number of pixel in X-direction")
    parser.add_option("--y_pix",dest="ypix", type=int, default = 600, help="number of pixel in Y-direction")
    parser.add_option("--rebin",dest="rebin", type=int, default = None, help="number of bins to group together")
    parser.add_option("--title",dest="iTitle", type=str, default = None, help="the title of the canvas")
    parser.add_option("--total",dest="total",action="store_true",default=False,help="use this flag to not plot individual source lines")
    parser.add_option("--catalogprefix",dest="catnames", type=str, default = "2FGL;P72Y", help="the prefix of the catalog sources, separate multiple names by ;")
    parser.add_option("--modelcounts",dest='modelcounts', action='store_true', default = False,
                      help = "use this flag to show integral over total model counts")
    (opts, arguments) = parser.parse_args()
    Input = glob.glob(sys.argv[1])
    iTitle = None
    Rootfile = None
    for i in Input:
        inputfile = str(i)
        if not opts.rootfile is None:
            Rootfile = opts.rootfile
            if len(Input)>1:
                Rootfile = inputfile.replace(".fits",".root")
        if not opts.iTitle is None:
            iTitle = opts.iTitle
        c = main(inputfile,rootfile=Rootfile,img_extension=opts.img_extension,xpix=opts.xpix,ypix=opts.ypix,iTitle=iTitle,catalog_names=opts.catnames,show_counts=opts.modelcounts,rebin=opts.rebin,plot_total=opts.total)
