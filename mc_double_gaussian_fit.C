using namespace RooFit;

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include "ACCSEL.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h> 
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooProduct.h>
#include <RooFormulaVar.h>
#include <RooExponential.h>
#include <RooGenericPdf.h>  
#include <RooProdPdf.h>

#include <RooFitResult.h>
#include <RooFit.h>
#include <RooCmdArg.h>



struct FitParams {
    double interval_fit_min, interval_fit_max;
    double mean_init, mean_min, mean_max;
    double sigma1_init, sigma1_min, sigma1_max;
    double sigma2_init, sigma2_min, sigma2_max;
    double c1_init, c1_min, c1_max;
    double N_min, N_max;
};


void fit_mc_signal_roofit(TString particle) {
    // Abrir ficheiro ROOT, obter a árvore, nome do eixo dos xx e definições do histograma, dependendo da partícula
    TTree* tree = nullptr;
    TString path_to_file = "";
    TString file_name = "";    
    TString path_to_tree = "";
    TString X_Axis_Title = "";
    TString legend_name = "";

    // Bins do Gráfico
    double xlow = 0.0;
    double xhigh = 0.0;
    double interval_of_bins = 0.0;
    int nbins = 0;

    // Posição das Legendas
    double x_min_legend = 0.0;
    double x_max_legend = 0.0;

    if (particle == "Bd"){
        path_to_file = "/lstore/cms/u25lekai/Bmeson/MC/ppRef/Bd_phat5_Bfinder.root";
        ///lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bd_phat5_Bfinder.root
        file_name = "Bd_Gaussian_Fit.pdf";
        path_to_tree = "Bfinder/ntKstar";
        X_Axis_Title = "m_{J/#Psi K^{*}} [GeV/c^{2}]";
        xlow = 5.0;
        xhigh = 6.0;
        interval_of_bins = 0.0075;
        nbins = int((xhigh - xlow) / interval_of_bins);
        x_min_legend = 0.50;
        x_max_legend = 0.88;
        legend_name = "B0";
        
    }
    else if (particle == "Bu"){
        path_to_file = "/lstore/cms/u25lekai/Bmeson/MC/ppRef/Bu_phat5_Bfinder.root";
        // /lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bu_phat5_Bfinder.root
        file_name = "Bu_Gaussian_Fit.pdf";
        path_to_tree = "Bfinder/ntKp";
        X_Axis_Title = "m_{J/#Psi K^{+}} [GeV/c^{2}]";
        xlow = 5.0;
        xhigh = 6.0;
        interval_of_bins = 0.0075;
        nbins = int((xhigh - xlow) / interval_of_bins);
        x_min_legend = 0.50;
        x_max_legend = 0.88;
        legend_name = "B^{+}";
    }
    else if (particle == "PSI2S"){
        path_to_file = "/lstore/cms/hlegoinha/X3872/MC_DATA/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder.root";
        file_name = "PSI2S_Gaussian_Fit.pdf";
        path_to_tree ="Bfinder/ntmix";
        X_Axis_Title = "m_{J/#Psi#pi^{+}#pi^{-}} [GeV/c^{2}]";
        xlow = 3.6;
        xhigh = 4.0;
        interval_of_bins = 0.0025;
        nbins = int((xhigh - xlow) / interval_of_bins);
        x_min_legend = 0.50;
        x_max_legend = 0.88;
        legend_name = "PSI2S";
    }
    else if (particle == "X3872"){
        path_to_file = "/lstore/cms/hlegoinha/X3872/MC_DATA/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder.root";
        file_name = "X3872_Gaussian_Fit.pdf";
        path_to_tree = "Bfinder/ntmix";
        X_Axis_Title = "m_{J/#Psi#rho} [GeV/c^{2}]";
        xlow = 3.6;
        xhigh = 4.0;
        interval_of_bins = 0.0025;
        nbins = int((xhigh - xlow) / interval_of_bins);
        x_min_legend = 0.14;
        x_max_legend = 0.52;
        legend_name = "X3872";
    }

    // Verificar se o Ficheiro e a Tree são encontrados
    TFile* file = TFile::Open(path_to_file);
    if (!file || file->IsZombie()) {
        std::cerr << "Ficheiro de MC não encontrado!" << std::endl;
        return;
    }

    file->GetObject(path_to_tree, tree);
    if (!tree) {
        std::cerr << "Tree não encontrada!" << std::endl;
        return;
    }

    // Definir cortes
    TString cut = "";
    if (particle == "Bu") {
        cut = Form("(%s) && (%s) && (%s) && (%s)",
                isMCsignal.Data(),
                ACCcuts_ppRef_Bu.Data(),
                SELcuts_ppRef_Bu.Data(),
                TRGmatching.Data());
    }
    else {
        cut = Form("(%s) && (%s) && (%s) && (%s)",
                isMCsignal.Data(),
                ACCcuts_ppRef.Data(),
                SELcuts_ppRef.Data(),
                TRGmatching.Data());
    }
    
    // Preencher histograma
    TH1F* h = new TH1F("h_roofit", "MC Signal;Bmass [GeV/c^{2}];Entries", nbins, xlow, xhigh);
    tree->Draw("Bmass >> h_roofit", cut, "goff");

    // Verificar se o histograma não fica vazio após os cortes (se fica, é sinal que algo está mal)
    if (h->GetEntries() == 0) {
        std::cerr << "Histograma vazio após os cortes!" << std::endl;
        delete h;
        return;
    }

    // Normalizar o histograma
    // h->Scale(1.0 / h->Integral());

    // VALORES INICIAIS E LIMITES DOS PARÂMETROS
    FitParams params;
    if (particle == "Bd"){
        params.interval_fit_min = 5.18;
        params.interval_fit_max = 5.38;

        params.mean_init = 5.28;
        params.mean_min = 5.26;
        params.mean_max = 5.3;

        params.sigma1_init = 0.014;
        params.sigma1_min = 0.001;
        params.sigma1_max = 0.07;

        params.sigma2_init = 0.034;
        params.sigma2_min = 0.001;
        params.sigma2_max = 0.09;

        params.c1_init = 0.8;
        params.c1_min  = 0.01;
        params.c1_max  = 0.99;
    
        params.N_min = 0;
        params.N_max = 90000;
    }
    else if (particle == "Bu"){
        params.interval_fit_min = 5.18;
        params.interval_fit_max = 5.42;

        params.mean_init = 5.28;
        params.mean_min = 5.26;
        params.mean_max = 5.3;

        params.sigma1_init = 0.036;
        params.sigma1_min = 0.001;
        params.sigma1_max = 0.09;

        params.sigma2_init = 0.011;
        params.sigma2_min = 0.001;
        params.sigma2_max = 0.07;

        params.c1_init = 0.52;
        params.c1_min  = 0.01;
        params.c1_max  = 0.99;   
    
        params.N_min = 0;
        params.N_max = 90000;
    }
    else if (particle == "PSI2S"){
        params.interval_fit_min = 3.64;
        params.interval_fit_max = 3.73;

        params.mean_init = 3.686;
        params.mean_min = 3.6;
        params.mean_max = 3.8;

        params.sigma1_init = 0.0048;
        params.sigma1_min = 0.001;
        params.sigma1_max = 0.09;

        params.sigma2_init = 0.0099;
        params.sigma2_min = 0.001;
        params.sigma2_max = 0.09;

        params.c1_init = 0.66;
        params.c1_min  = 0.1;
        params.c1_max  = 0.9;
    
        params.N_min = 0;
        params.N_max = 100000;
    }
    else if (particle == "X3872"){
        params.interval_fit_min = 3.8;
        params.interval_fit_max = 3.95;

        params.mean_init = 3.872;
        params.mean_min = 3.8;
        params.mean_max = 3.9;

        params.sigma1_init = 0.0056;
        params.sigma1_min = 0.001;
        params.sigma1_max = 0.08;

        params.sigma2_init = 0.0125;
        params.sigma2_min = 0.001;
        params.sigma2_max = 0.08;

        params.c1_init = 0.473;
        params.c1_min  = 0.01;
        params.c1_max  = 0.99;
    
        params.N_min = 0;
        params.N_max = 100000;
    }

    // RooFit: variável observável
    RooRealVar Bmass("Bmass", "Bmass", xlow, xhigh);
    Bmass.setRange("fitRange", params.interval_fit_min, params.interval_fit_max);
    RooDataHist data("data", "Histograma RooFit", RooArgList(Bmass), h);

    // Parâmetros do fit
    RooRealVar mean("mean", "mean", params.mean_init, params.mean_min, params.mean_max);
    RooRealVar sigma1("sigma1", "sigma1", params.sigma1_init, params.sigma1_min, params.sigma1_max);
    RooRealVar sigma2("sigma2", "sigma2", params.sigma2_init, params.sigma2_min, params.sigma2_max);
    RooRealVar c1("c1", "Fraction of G1", 0.65, 0.0, 1.0);

    RooRealVar N("N", "Total yield", h->GetEntries(), params.N_min, params.N_max);

    // Gaussianas
    RooGaussian gauss1("gauss1", "Narrow", Bmass, mean, sigma1);
    RooGaussian gauss2("gauss2", "Wide", Bmass, mean, sigma2);

    // Soma normalizada
    RooAddPdf shape("shape", "Normalized Double Gaussian", RooArgList(gauss1, gauss2), RooArgList(c1));

    // Modelo estendido com yield total
    RooExtendPdf model("model", "Extended model", shape, N);

    // Fit
    RooFitResult* result = model.fitTo(data, Range("fitRange"), Save(), Extended(kTRUE));

    // Valores Legenda
    double mean_val = mean.getVal();
    double mean_err = mean.getError();
    double N_val = N.getVal();
    double N_err = N.getError();

    //std::cout << "CANDIDATES: " << h->GetEntries() << std:endl;

    // Desenhar
    TCanvas* c = new TCanvas("c_roofit", "RooFit", 800, 600);
    RooPlot* frame = Bmass.frame();
    frame->SetTitle(" ");

    // Dar Nome ao Eixo dos xx Dependendo da Partícula
    frame->GetXaxis()->SetTitle(X_Axis_Title);
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", interval_of_bins));
    frame->GetYaxis()->SetTitleOffset(1.35);

    data.plotOn(frame, MarkerStyle(20), MarkerSize(1.0), Name("h_roofit")); // Circle marker
    model.plotOn(frame, LineColor(kBlue), Name("model"), Range("fitRange"), NormRange("fitRange")); // Blue model line
    shape.plotOn(frame, Components(gauss1), LineStyle(kDashed), LineColor(kGreen + 2), Range("fitRange"), NormRange("fitRange"));
    shape.plotOn(frame, Components(gauss2), LineStyle(kDashed), LineColor(kRed), Range("fitRange"), NormRange("fitRange"));
    
    frame->Draw();

    int nFitParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("model", "h_roofit", nFitParams);
    
    // Legend 1 — MC Signal, Model, Number of Candidates
    auto legend_signal = new TLegend(x_min_legend, 0.76, x_max_legend, 0.88);
    legend_signal->SetTextSize(0.03);
    legend_signal->SetTextFont(42);
    legend_signal->AddEntry((TObject*)frame->findObject("h_roofit"), legend_name + " - MC Signal", "p");
    legend_signal->AddEntry((TObject*)frame->findObject("model"), "Model: Double Gaussian", "l");
    legend_signal->AddEntry((TObject*)0, Form("No. of MC Candidates = %.0f", h->GetEntries()), "");
    legend_signal->Draw();

    // TPaveText for fit parameters (more compact, no indent)
    TPaveText* pave = new TPaveText(x_min_legend, 0.47, x_max_legend, 0.75, "NDC");
    pave->SetTextAlign(12); // left align
    pave->SetTextFont(42);
    pave->SetTextSize(0.03);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);

    pave->AddText(Form("Mean = %.5f #pm %.5f", mean.getVal(), mean.getError()));
    pave->AddText(Form("#sigma_{1} = %.5f #pm %.5f", sigma1.getVal(), sigma1.getError()));
    pave->AddText(Form("#sigma_{2} = %.5f #pm %.5f", sigma2.getVal(), sigma2.getError()));
    pave->AddText(Form("c_{1} = %.4f #pm %.4f", c1.getVal(), c1.getError()));
    pave->AddText(Form("N = %.2f #pm %.2f", N.getVal(), N.getError()));
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
    pave->Draw();

    // Terminal
    std::cout << " ➤ Mean = " << mean.getVal() << std::endl;
    std::cout << " ➤ Sigma1 = " << sigma1.getVal() << std::endl;
    std::cout << " ➤ Sigma2 = " << sigma2.getVal() << std::endl;
    std::cout << " ➤ c1 (G1 fraction) = " << c1.getVal() << std::endl;
    std::cout << " ➤ N = " << N_val << std::endl;
    std::cout << " ➤ Total entries = " << h->GetEntries() << std::endl;

    // Guardar com Nome Específico da Partícula
    c->SaveAs(file_name);

    // Limpar
    delete h;
    delete c;
}


// Particles 2G Monte Carlo: Bd, Bu, X3872, PSI2S
void mc_double_gaussian_fit() {
    fit_mc_signal_roofit("Bd");
    //fit_mc_signal_roofit("Bu");
    //fit_mc_signal_roofit("X3872");
    //fit_mc_signal_roofit("PSI2S");
}