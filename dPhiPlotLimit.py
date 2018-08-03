#!/usr/bin/env python

import os, glob, sys, math
from array import array
from ROOT import *
import numpy as np
#def GetContours(g, color, style):
#    contours = [1.0]
#    g.GetHistogram().SetContour(1,array('d',contours));
#    g.Draw("cont z list");
#    contLevel = g.GetContourList(1.0);
#    max_points = -1
#    for cont in contLevel:
#        n_points = cont.GetN()
#        if n_points > max_points:
#            max_points = n_points
#            cont.SetLineColor(color)
#            cont.SetLineStyle(style)
#            cont.SetLineWidth(5)
#            out = cont
#
#            #cont.Draw("same")
#
#    return out



def GetContoursSmooth(g, color, style, n_smooth  = 4, limit = True):

    if(n_smooth>0):
        g2 = g.Clone()
        histo2d = g2.GetHistogram();
        htemp = TH2D("", "",
                100, histo2d.GetXaxis().GetXmin(), histo2d.GetXaxis().GetXmax(),
                100, histo2d.GetYaxis().GetXmin(), histo2d.GetYaxis().GetXmax());
        for binx in range(1,htemp.GetNbinsX()):
            x = htemp.GetXaxis().GetBinCenter(binx)
            for biny in range(1,htemp.GetNbinsY()):
                y = htemp.GetYaxis().GetBinCenter(biny)
                z = g2.Interpolate(x,y);
                if(z!=0.):
                    htemp.SetBinContent(htemp.GetBin(binx, biny), z)

        for ind in range(0,n_smooth):
            htemp.Smooth(1,"k5b");


        vx=[]; vy=[]; vz=[];
        glu_lsp = 225
        for binx in range(1,htemp.GetNbinsX()):
            x = htemp.GetXaxis().GetBinCenter(binx)
            for biny in range(1,htemp.GetNbinsY()):
                y = htemp.GetYaxis().GetBinCenter(biny);
                z = htemp.GetBinContent(htemp.GetBin(binx,biny));

                vx.append(x)
                vy.append(y)
                if ((x-y) > (glu_lsp+85)):
                    vz.append(z)
                else:
                    vz.append(g2.Interpolate(x,y));
        ax = array("d", vx)
        ay = array("d", vy)
        az = array("d", vz)
        gsmooth =  TGraph2D ("gsmooth", "Cross-Section Limit", len(vx), ax, ay, az)
    else:
        gsmooth = g.Clone()

    if limit:
        contours = [1.0]
    else:
        print 'Foo08'
        contours = [5.]
    #contours = [1.,2.,3.,4.,5.]
    gsmooth.GetHistogram().SetContour(1,array('d',contours));
    gsmooth.Draw("cont z list");
    contLevel = gsmooth.GetContourList(contours[0]);
    max_points = -1
    #find the contour with the most points
    outSM = TGraph()
    for i,cont in enumerate(contLevel):
        n_points = cont.GetN()
        if n_points > max_points:
            max_points = n_points
            outSM = cont

    #get the unsmoothed contour anyway for the last diagonal point
    #contours = [1.0]
    g.GetHistogram().SetContour(1,array('d',contours));
    contLevel = g.GetContourList(contours[0]);
    max_points = -1
    outnSM = TGraph()
    #find the contour with the most points
    for i,cont in enumerate(contLevel):
        n_points = cont.GetN()
        if n_points > max_points:
            max_points = n_points
            outnSM = cont


    outnSM.SetLineColor(color)
    outnSM.SetLineStyle(style)
    outnSM.SetLineWidth(5)

    glu_lsp = 225
    #First: remove line above the diagonal
    if n_smooth > 0:
        for point in range(0,outSM.GetN()):
            mglu, mlsp = Double(0), Double(0)
            outSM.GetPoint(point, mglu, mlsp);
            if(mlsp > mglu-glu_lsp-5):
                while(point <= outSM.GetN() and point!=0):
                    outSM.RemovePoint(outSM.GetN()-1)
    #Second:extend line down to LSP =0
#    outSM.Sort() #sorting doesn't really work
        endglu, endlsp = Double(0), Double(0)
        outSM.GetPoint(1, endglu, endlsp)
        outSM.SetPoint(1, endglu, 0)

    #Extend line on the diagonal
#        iniglu, inilsp = Double(0), Double(0)
       # outSM.GetPoint(outSM.GetN()-1, iniglu, inilsp)
       # outSM.SetPoint(outSM.GetN()-1, inilsp+225, inilsp)



    outSM.SetLineColor(color)
    outSM.SetLineStyle(style)
    outSM.SetLineWidth(5)
    return outSM

def getxsecGlu():
    xsecGlu = {} # dict for xsecs
    #xsecFile = "../../../../../SUSYAnalysis/python/tools/glu_xsecs_13TeV.txt"
    xsecFile = "./tchiwz_xs_14tev.txt"

    with open(xsecFile,"r") as xfile:
        lines = xfile.readlines()
        print 'Found %i lines in %s' %(len(lines),xsecFile)
        for line in lines:
            if line[0] == '#': continue
            (mGo,xsec,err) = line.split()
            xsecGlu[int(mGo)] = (float(xsec),float(err))
            #print 'Importet', mGo, xsec, err, 'from', line
    return xsecGlu

if __name__ == "__main__":

    ## remove '-b' option
    if '-b' in sys.argv:
        sys.argv.remove('-b')
        _batchMode = True

    if len(sys.argv) > 1:
        pattern = sys.argv[1]
        print '# pattern is', pattern
    else:
        print "No pattern given!"
        exit(0)

    ## Create Yield Storage

#    pattern = "datacardsABCD_2p1bins_fullscan2"
#    os.chdir(pattern)
    #dirList = glob.glob('T1tttt*')
    #samples = [x[x.find('/')+1:] for x in dirList]
    xsecGlu = getxsecGlu()
    if sys.argv[2] == "limit":
        limit = True
    else:
        limit = False # if False, significance
    if 1==1:
        #fileList = glob.glob(pattern+'/limitOutput/*root')
        if limit:
            fileList = glob.glob('../HiggsAnalysis/CombinedLimit/scripts/limits/limit_*.root')
        else:
            fileList = glob.glob('../HiggsAnalysis/CombinedLimit/scripts/limits/significance_*.root')
        #hexp = TH2F('hexp','hexp', 81,-12.5, 2012.5, 81, -12.5, 2012.5)
        #hexpdown = TH2F('hexpdown','hexpdown', 81,-12.5, 2012.5, 81, -12.5, 2012.5)
        #hexpup = TH2F('hexpup','hexpup', 81,-12.5, 2012.5, 81, -12.5, 2012.5)
        #hobs = TH2F('hobs','hobs', 81,-12.5, 2012.5, 81, -12.5, 2012.5)
        hexp = TH2F('hexp','hexp', 11, -25., 525., 11, -2.5, 52.5)
        hexpdown = TH2F('hexpdown','hexpdown', 11, -25., 525., 11, -2.5, 52.5)
        hexpup = TH2F('hexpup','hexpup', 11, -25., 525., 11, -2.5, 52.5)
        hobs = TH2F('hobs','hobs', 11, -25., 525., 11, -2.5, 52.5)
        vmx=[]; vmy = []; vxsec = []; vobs = [];  vobsup = []; vobsdown = []; vexp = []; vup = []; vdown = []; vlim = [];
        for x in fileList:
            mGo = int(x.split('_')[2])
            mLSP = int(x.split('_')[3].rstrip('.root'))
            #mGo = int(x[x.find('_mGo')+4:x.find('_mLSP')])
            #mLSP = int(x[x.find('_mLSP')+5:x.find('.As')])
            f = TFile.Open(x, 'read')
            t = f.Get('limit')
            xsec = xsecGlu[mGo][0]
            theorySys = xsecGlu[mGo][1]
            #print theorySys
            rExp = 0
            rObs = 0
            rExp1SigmaDown = 0
            rExp1SigmaUp = 0
            factor = 1.0
            #if mGo < 1400:
            #    factor = 100.0
            if not t:
                print 'Cannot find tree in file {}. Skip.'.format(x)
                continue
            for entry in t:
                q = entry.quantileExpected
                if q == 0.5: rExp = entry.limit/factor
                if q == -1: rObs = entry.limit/factor
                if q < 0.4 and q > 0.14 : rExp1SigmaDown = entry.limit/factor
                if q < 0.14 : rExp2SigmaDown = entry.limit/factor
                if q > 0.6 and q < 0.9 : rExp1SigmaUp = entry.limit/factor
                if q > 0.9 : rExp2SigmaUp = entry.limit/factor
            if limit == False:
                rExp = rObs
            print 'Foo01', mGo, mLSP, rExp, rObs

            hexp.Fill(mGo,mLSP,rExp)
            hexpdown.Fill(mGo,mLSP,rExp1SigmaDown)
            hexpup.Fill(mGo,mLSP,rExp1SigmaUp)
            hobs.Fill(mGo,mLSP,rObs)
            vmx.append(mGo)
            vmy.append(mLSP)
            vxsec.append(xsec)
            vlim.append(xsec * rObs)
            vobs.append(rObs)
            vobsup.append(rObs*(1+theorySys/100.0))
            vobsdown.append(rObs*(1-theorySys/100.0))
            vexp.append(rExp)
            vup.append(rExp1SigmaUp)
            vdown.append(rExp1SigmaDown)
            f.Close()

        hexp.SaveAs(pattern+'/testexp_'+pattern+'.root')
        hobs.SaveAs(pattern+'/testobs_'+pattern+'.root')
        aexp = array("d", vexp)
        alim = array("d", vlim)
        aup = array("d", vup)
        adown = array("d", vdown)
        aobs = array("d", vobs)
        aobsup = array("d", vobsup)
        aobsdown = array("d", vobsdown)
        amx = array("d", vmx)
        amy = array("d", vmy)

        glim = TGraph2D("glim", "Cross-section limt", len(vlim), amx, amy, alim)
        gexp = TGraph2D("gexp", "Expected Limit", len(vexp), amx, amy, aexp)
        gup = TGraph2D("gup", "Expected Limit 1sigma up", len(vup), amx, amy, aup)
        gdown = TGraph2D("gdown", "Expected Limit 1sigma down", len(vdown), amx, amy, adown)
        gobs = TGraph2D("gobs", "Observed Limit", len(vobs), amx, amy, aobs)
        gobsup = TGraph2D("gobsup", "theory 1sigma up", len(vobsup), amx, amy, aobsup)
        gobsdown = TGraph2D("gobsdown", "theory 1sigma down", len(vobsdown), amx, amy, aobsdown)
        c = TCanvas("c","c",800,600)

        xmin = min(vmx)
        xmax = max(vmx)
        ymin = min(vmy)
        ymax = max(vmy)
        bin_size = 12.5;
        nxbins = max(1, min(500, (math.ceil((xmax-xmin)/bin_size))))
        nybins = max(1, min(500, (math.ceil((ymax-ymin)/bin_size))))
        glim.SetNpx(int(nxbins))
        glim.SetNpy(int(nybins))


        cexp = GetContoursSmooth(gexp, 2,1,0, limit)
        #cup = GetContoursSmooth(gup,2,2,0)
        #cdown = GetContoursSmooth(gdown,2,2,0)
        #cobs = GetContoursSmooth(gobs,1,1,0)
        #cobsup = GetContoursSmooth(gobsup,1,2,0)
        #cobsdown = GetContoursSmooth(gobsdown,1,2,0)

        #cexp = GetContoursSmooth(gexp, 2,1)
        #cup = GetContoursSmooth(gup,2,2)
        #cdown = GetContoursSmooth(gdown,2,2)
        #cobs = GetContoursSmooth(gobs,1,1)
        #cobsup = GetContoursSmooth(gobsup,1,2)
        #cobsdown = GetContoursSmooth(gobsdown,1,2)

        hlim = glim.GetHistogram()
        hlim.SetTitle(";m_{gluino} [GeV];m_{LSP} [GeV]");
        hlim.Draw("colz")
        cexp.Draw("same")
        x = np.ndarray((cexp.GetN()), 'd', cexp.GetX())
        y = np.ndarray((cexp.GetN()), 'd', cexp.GetY())
        print 'Foo05', list(x)
        print 'Foo06', list(y)
        ### FOO HACK
        listx = list(x)
        listy = list(y)
        #print 'Foo07', len(listx), len(listy)
        #for idx in range(0, len(listx)):
        #    print int(listx[idx]), int(listy[idx])
        #print 'Foo08'
        ## SMOOTH
        #listx[1] -= 10.
        #listx[2] += 10.
        #listx[12] = listx[11]
        #listx[13] = listx[11]
        #listx[14] = listx[11]
        #listy[12] = listy[11]
        #listy[13] = listy[11]
        #listy[14] = listy[11]
        #listx.append(196)
        #listy.append(43)
        #print 'Foo09', len(listx), len(listy)
        #for idx in range(0, len(listx)):
        #    print int(listx[idx]), int(listy[idx])
        #print 'Foo10'
        ##listx.append(250.)
        ##listy.append(50.)
        newx = array("d", listx)
        newy = array("d", listy)
        #newx = array("d", [80., 110., 120., 144., 158.0, 153.6, 153.4, 151.4, 150.4, 132.0, 110.0, 100.])
        #newy = array("d", [5., 8., 9., 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 39.0, 40.0, 42.3])
        cexp = TGraph(len(newx), newx, newy)

        #cup.Draw("same")
        #cdown.Draw("same")
        #cobs.Draw("same")
        #cobsup.Draw("same")
        #cobsdown.Draw("same")
        flimit = TFile(pattern+"/limit_scan.root","recreate")

        #cobs.SetTitle("Observed Limit");
        #cobsup.SetTitle("Observed -1#sigma Limit");
        #cobsdown.SetTitle("Observed +1#sigma Limit");
        cexp.SetTitle("Expected  Limit");
        #cup.SetTitle("Expected -1#sigma Limit");
        #cdown.SetTitle("Expected +1#sigma Limit");


        hlim.Write("T1ttttObservedExcludedXsec");
        #cobs.Write("T1ttttObservedLimit");
        #cobsup.Write("T1ttttObservedLimitDown");
        #cobsdown.Write("T1ttttObservedLimitUp");
        cexp.Write("T1ttttExpectedLimit");
        #cup.Write("T1ttttExpectedLimitDown");
        #cdown.Write("T1ttttExpectedLimitUp");

        c.SaveAs(pattern+'/canvas_'+pattern+'.root')
        c.SaveAs(pattern+'/canvas_'+pattern+'.pdf')

