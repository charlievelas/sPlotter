#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooStats/SPlot.h>
#include <RooPolynomial.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLegend.h>
#include <iostream>
#include <chrono>

void sPlot_dev3() {

    // Load RooFit gSystem
    gSystem->Load("libRooFit");
    gSystem->Load("libRooStats.so");

    // Load root file and tree
    TFile f_1("distributions_order2.root");
    TTree *tree_1 = (TTree*)f_1.Get(f_1.GetListOfKeys()->At(0)->GetName());

    // Define the names of the branches
    std::string x_name = "mass";
    std::string y_name = "pz";
    int pzAIndx;
    tree_1->SetBranchAddress("pzAIndx",&pzAIndx);

    // Define the RooRealVar variables corresponding to the TTree branches
    RooRealVar x(x_name.c_str(), x_name.c_str(), 0, 10); // Assuming the branch "mass" has values between 0 and 10
    RooRealVar y(y_name.c_str(), y_name.c_str(), -6, 6); // Assuming the branch "pz" has values between -5 and 5

    // Create a RooArgSet to hold the variables
    RooArgSet vars(x, y);

    // Create the RooDataSet from the TTree
    RooDataSet* data = new RooDataSet("data", "dataset from TTree", tree_1, vars);

    // Print data
    std::cout << "data printed " << std::endl;
    data->Print("v");

    // Create a RooDataSet containing only the x data
    RooDataSet* data_x = (RooDataSet*)data->reduce(RooArgSet(x));

    // Create signal and background yields - initially assume half signal and half background
    int num_of_evnts = data_x->numEntries();
    RooRealVar nsig("nsig", "number of signal events", num_of_evnts/2, 0, num_of_evnts);
    RooRealVar nbkg("nbkg", "number of background events", num_of_evnts/2, 0, num_of_evnts);

    // Define discriminating variable PDFs and combine
    RooRealVar mean("mean", "mean of gaussian", 5, 0, 10);
    RooRealVar sigma("sigma", "width of gaussian", 1, 0.1, 10);
    RooRealVar a0("a0", "a0", 0.1,-1,1);
    RooRealVar a1("a1", "a1", 0.01,-1,1);
    RooRealVar a2("a2", "a2", 0.01,-1,1);
    RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);
    RooChebychev bkg("bkg", "background PDF", x, RooArgList(a0,a1,a2));
    RooAddPdf model_x("model_x", "signal + background", RooArgList(gauss, bkg), RooArgList(nsig, nbkg));

    // Fit model to data_x
    RooFitResult *fitResult = model_x.fitTo(*data_x, RooFit::Save());

    // Print covariance matrix
    TMatrixDSym cov_matrix = fitResult->covarianceMatrix();
    std::cout << "Printing covariance matrix after fitting:" << std::endl;
    cov_matrix.Print();

    // Print the signal and background yields
    std::cout << "Signal yield (x): " << nsig.getVal() << " ± " << nsig.getError() << std::endl;
    std::cout << "Background yield (x): " << nbkg.getVal() << " ± " << nbkg.getError() << std::endl;

    // Create sPlot object for data_x
    RooStats::SPlot* sData = new RooStats::SPlot("sData", "An SPlot", *data, &model_x, RooArgList(nsig, nbkg));

    // Create sWeighted dataset for y
    RooDataSet* dataw_sig = new RooDataSet("dataw_sig", "sWeighted data_y", data, *data->get(), 0, "nsig_sw");
    RooDataSet* dataw_bg = new RooDataSet("dataw_bg", "sWeighted data_y", data, *data->get(), 0, "nbkg_sw");

    // Create canvas with 4 pads
    TCanvas* c = new TCanvas("sPlot input", "sPlot input", 1200, 1200);
    c->Divide(2, 2);

    // Pad 1: Original distribution with fit
    c->cd(1);
    RooPlot* frame1 = x.frame();
    data->plotOn(frame1);
    model_x.plotOn(frame1);
    model_x.plotOn(frame1, RooFit::Components(bkg), RooFit::LineStyle(kDashed));
    model_x.plotOn(frame1, RooFit::Components(gauss), RooFit::LineColor(kRed));
    frame1->SetTitle("Original Distribution with Fit");
    frame1->GetYaxis()->SetTitleOffset(1.2); // Move y-axis label closer to the axis
    frame1->Draw();

    // Pad 2: Signal and Background sWeights vs x
    c->cd(2);
    TGraph* graph_sWeights_sig = new TGraph();
    TGraph* graph_sWeights_bkg = new TGraph();
    TGraph* graph_sWeights_sum = new TGraph();
    double minWeight = 0;
    double maxWeight = 0;
    std::vector<double> s_weights;
    std::vector<double> b_weights;
    for (int i = 0; i < data_x->numEntries(); ++i) {
        const RooArgSet* row = data_x->get(i);
        double x_val = row->getRealValue(x_name.c_str());
        double sWeight_sig = sData->GetSWeight(i, "nsig");
        double sWeight_bkg = sData->GetSWeight(i, "nbkg");
        double sWeight_sum = sWeight_sig + sWeight_bkg;
        graph_sWeights_sig->SetPoint(i, x_val, sWeight_sig);
        graph_sWeights_bkg->SetPoint(i, x_val, sWeight_bkg);
        graph_sWeights_sum->SetPoint(i, x_val, sWeight_sum);
        if (sWeight_sig < minWeight) minWeight = sWeight_sig;
        if (sWeight_bkg < minWeight) minWeight = sWeight_bkg;
        if (sWeight_sum < minWeight) minWeight = sWeight_sum;
        if (sWeight_sig > maxWeight) maxWeight = sWeight_sig;
        if (sWeight_bkg > maxWeight) maxWeight = sWeight_bkg;
        if (sWeight_sum > maxWeight) maxWeight = sWeight_sum;
        s_weights.push_back(sWeight_sig);
        b_weights.push_back(sWeight_bkg);
    }
    graph_sWeights_sig->SetTitle("Signal and Background sWeights vs x; x; sWeight");
    graph_sWeights_sig->SetMarkerStyle(1);
    graph_sWeights_sig->SetMarkerColor(kRed);
    graph_sWeights_sig->GetXaxis()->SetLimits(0, 10); // Ensure x-axis range matches mass x-axis range
    graph_sWeights_sig->GetXaxis()->SetTitle(x_name.c_str());
    graph_sWeights_sig->GetYaxis()->SetTitle("sWeight");
    graph_sWeights_sig->GetYaxis()->SetTitleOffset(1.2); // Move y-axis label closer to the axis
    graph_sWeights_sig->Draw("AP");
    graph_sWeights_bkg->SetMarkerStyle(1);
    graph_sWeights_bkg->SetMarkerColor(kBlue);
    graph_sWeights_bkg->Draw("P SAME");
    graph_sWeights_sum->SetMarkerStyle(1);
    graph_sWeights_sum->SetMarkerColor(kBlack);
    graph_sWeights_sum->Draw("P SAME");

    // Set y-axis range to ensure all weights are plotted
    graph_sWeights_sig->GetYaxis()->SetRangeUser(minWeight - 0.1 * fabs(minWeight), maxWeight + 0.1 * fabs(maxWeight));

    // Pad 3: Plot data_y
    c->cd(3);
    RooPlot* frame3 = y.frame();
    data->plotOn(frame3);
    frame3->SetTitle("Distribution of y");
    frame3->GetYaxis()->SetTitleOffset(1.2); // Move y-axis label closer to the axis
    frame3->Draw();
    tree_1->Draw("pz>>h_s(100,-6,6)", "pzAIndx==1", "SAME");
    TH1F *h_s = (TH1F*)gDirectory->Get("h_s");
    h_s->SetLineColor(kBlue);
    h_s->SetLineStyle(2);
    h_s->Draw("SAME");
    tree_1->Draw("pz>>h_b(100,-6,6)", "pzBIndx==1", "SAME");
    TH1F *h_b = (TH1F*)gDirectory->Get("h_b");
    h_b->SetLineColor(kRed);
    h_b->SetLineStyle(2);
    h_b->Draw("SAME");

    // Pad 4: Plot sWeighted data_y
    c->cd(4);
    RooPlot* frame_sig_y = y.frame();
    data->plotOn(frame_sig_y);
    // Draw the original pz distribution when pzAIndx == 1 on the same pad
    // Plot sWeighted pz
    dataw_sig->plotOn(frame_sig_y, RooFit::DataError(RooAbsData::SumW2), RooFit::LineColor(kRed));
    frame_sig_y->SetTitle("sWeighted Distribution of y");
    frame_sig_y->GetYaxis()->SetTitleOffset(1.2); // Move y-axis label closer to the axis
    frame_sig_y->Draw();
    h_s->Draw("SAME");

    // Save the canvas
    c->Print("sPlot_dev3.root");

    // Create canvas for cumulative probability distributions
    TCanvas* cr = new TCanvas("sPlot ratio", "sPlot ratio", 1200, 600);
    TH1F *h_cdf_s = (TH1F*)h_s->GetCumulative();
    TH1F *h_cdf_w = new TH1F("h_cdf_w", "Cumulative Distribution of sWeighted pz", 100, -6, 6);
    for (int i = 0; i < dataw_sig->numEntries(); ++i) {
        const RooArgSet* row = dataw_sig->get(i);
        double pz_val = row->getRealValue("pz");
        double weight = dataw_sig->weight();
        h_cdf_w->Fill(pz_val, weight);
    }
    h_cdf_w = (TH1F*)h_cdf_w->GetCumulative();

    h_cdf_s->SetLineColor(kBlue);
    h_cdf_s->SetLineWidth(2);
    h_cdf_s->SetTitle("Cumulative Probability Distributions");
    h_cdf_s->Draw();

    h_cdf_w->SetLineColor(kRed);
    h_cdf_w->SetLineWidth(2);
    h_cdf_w->Draw("SAME");

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h_cdf_s, "True Distribution", "l");
    legend->AddEntry(h_cdf_w, "sWeighted Distribution", "l");
    legend->Draw();

    // Perform Kolmogorov-Smirnov test
    double ksTest = h_s->KolmogorovTest(h_cdf_w);
    std::cout << "Kolmogorov-Smirnov Test p-value = " << ksTest << std::endl;

    // Save the second canvas
    cr->Print("sPlot_dev3_ratio.root");

    // Save sWeights to a new ROOT file
    std::string input_filename = f_1.GetName();
    std::string output_filename = input_filename.substr(0, input_filename.find_last_of('.')) + "_sWeights.root";
    TFile output_file(output_filename.c_str(), "RECREATE");
    // Clone the original input tree
    TTree* cloned_tree = tree_1->CloneTree(0);
    // Create new branches for signal and background weights
    float sWeights_s;
    float sWeights_b;
    cloned_tree->Branch("sWeights_s", &sWeights_s, "sWeights_s/F");
    cloned_tree->Branch("sWeights_b", &sWeights_b, "sWeights_b/F");
    // Fill the cloned tree with the new branches
    for (size_t i = 0; i < s_weights.size(); ++i) {
        tree_1->GetEntry(i);
        sWeights_s = s_weights[i];
        sWeights_b = b_weights[i];
        cloned_tree->Fill();
    }
    // Save the cloned tree with the new branches to the output file
    output_file.cd();
    cloned_tree->Write();
    output_file.Close();
}

int main() {
    sPlot_dev3();
    return 0;
}