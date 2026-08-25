void Fit() {
	// 1. Clone histogram in memory so original remains intact
	TH2F *h2d_gated = (TH2F*)p_tail_shape->Clone("p_tail_shape_gated"); // For protons
	//TH2F *h2d_gated = (TH2F*)a_tail_shape->Clone("a_tail_shape_gated"); // For alphas

	// 2. Linear cut points

	// For proton fit
	const Double_t x1 = -168.;
	const Double_t y1 = 680.;
	const Double_t x2 = 0.;
	const Double_t y2 = 718.;

	// For alpha fit
	//const Double_t x1 = -177.;
	//const Double_t y1 = 492.;
	//const Double_t x2 = 0.;
	//const Double_t y2 = 530.;

	// 3. Zero out gated bins below the line
	for (int i = 1; i <= h2d_gated->GetNbinsX(); ++i) {
		Double_t x_val = h2d_gated->GetXaxis()->GetBinCenter(i);
		Double_t y_threshold = y1 + ((y2 - y1) / (x2 - x1)) * (x_val - x1);

		for (int j = 1; j <= h2d_gated->GetNbinsY(); ++j) {
			Double_t y_val = h2d_gated->GetYaxis()->GetBinCenter(j);
			if (y_val < y_threshold) {
				h2d_gated->SetBinContent(i, j, 0);
				h2d_gated->SetBinError(i, j, 0);
			}
		}
	}

	TF1 *fit1D = new TF1("fit1D", "[0]*exp([1]*x)");
	fit1D->SetParameters(567.0, 0.00034); // initial guesses from non-gated 1D fit

	h2d_gated->Fit(fit1D);
}
