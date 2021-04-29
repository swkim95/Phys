#include "./setTDRStyle.C"

void getControlPlots( void )
{
  TString _channel = "mumu";
  TString _opt = "mass";
  TString _ptCut = "_pt30";

  gStyle->SetOptStat(0);

  TFile* f1 = new TFile("./rootfiles/massSpec_"+_channel+"_"+_opt+_ptCut+".root");
  f1->cd();
  // data
  TH1D* hdata;
  hdata = (TH1D*) hist_SingleMu_Run2011B->Clone();
  // MC
  TH1D* hmc;
  hmc = (TH1D*) hist_DYJets->Clone();
  cout << "nevents = (data): " << hdata->Integral() << ";; (DY MC): " << hmc->Integral() << endl;
  // Normalize to the number of data
  hmc->Scale(hdata->Integral()/hmc->Integral());

  hmc->SetLineColor(kYellow);
  hmc->SetFillColor(kYellow);
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->cd();
  c1->SetLogy();
  TPad *c1_1 = new TPad("padc1_1","padc1_1",0.01,0.05,0.99,0.99);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetTopMargin(0.01);
  c1_1->SetBottomMargin(0.3);
  c1_1->SetRightMargin(0.1);
  c1_1->SetLeftMargin(0.109);
  c1_1->SetFillStyle(0);
  //c1_1->SetLogy();
  
  hdata->SetTitle("CMS Preliminary at #sqrt{s} = 7 TeV, 2011 data");
  hdata->Draw("ep");
  hdata->SetLabelSize(0.0);
  hdata->GetYaxis()->SetTitleOffset(1.5);
  hdata->GetYaxis()->SetTitle("Events / 1 GeV");
  hdata->SetMarkerStyle(20);
  hmc->Draw("histsame");
  hdata->Draw("epsame");

  ostringstream tmpstrm;
  TLegend *legend = new TLegend(0.70, 0.70, 0.90, 0.90);
  tmpstrm.str("");
  tmpstrm << "DATA   " ;
  legend->AddEntry(hdata, tmpstrm.str().c_str(), "PL");
  tmpstrm.str("");
  tmpstrm << "#gamma^{*}/Z #rightarrow #mu#mu   " ;
  legend->AddEntry(hmc, tmpstrm.str().c_str(), "F");
  tmpstrm.str("");

  legend->SetFillColor(0);
  legend->Draw("0");

  c1_2 = new TPad("padc1_2","padc1_2",0.01,0.05,0.99,0.32);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetTopMargin(0.1);
  c1_2->SetBottomMargin(0.30);
  c1_2->SetRightMargin(0.091);
  //c1_2->SetLeftMargin(0.122);
  c1_2->SetFillStyle(0);
  c1_2->SetGrid();
  //c1_2->SetLogx();
  TH1D* hratio = (TH1D*) hdata->Clone();
  hdata->Sumw2(); hmc->Sumw2();
  hratio->Divide(hdata, hmc);
  hratio->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
  hratio->SetTitle("");
  hratio->GetXaxis()->SetMoreLogLabels();
  hratio->GetXaxis()->SetNoExponent();
  hratio->GetYaxis()->SetTitle("data/MC");
  hratio->GetXaxis()->SetTitleSize(0.13);
  hratio->GetYaxis()->SetTitleSize(0.09);
  hratio->GetYaxis()->SetTitleOffset(0.4);
  hratio->GetXaxis()->SetLabelSize(0.11);
  hratio->GetYaxis()->SetLabelSize(0.07);
  
  gStyle->SetOptFit(0);
  hratio->SetMaximum(2.0);
  hratio->SetMinimum(0.0);
  hratio->SetMarkerSize(0.5);
  hratio->Draw("e1p");
  hratio->Fit("pol0");

  //c1->Print("./eps/plot_"+_channel+"_"+_opt+".png");
}
