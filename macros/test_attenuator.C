TGraph*       NewGraph      (TString name="graph", int mst=20, double msz=0.6, int mcl=kBlack, int lst=-1, int lsz=-1, int lcl=-1);
TGraphErrors* NewGraphErrors(TString name="graph", int mst=20, double msz=0.6, int mcl=kBlack, int lst=-1, int lsz=-1, int lcl=-1);
double CircleIntersectionArea(double x1, double y1, double r1, double x2, double y2, double r2);
double Gaussian2D(double x, double y, double xc, double yc, double sigmaX, double sigmaY);
double IntegrateGaussian2D(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints=100);
double IntegrateGaussian2D_SumGamma(double sigma, double cx, double r0);
double ApproximateGaussianIntegral(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0);

void test_attenuator()
{
    double wFull = 40; // width of active tantalium aread
    double xMid = 0.5*wFull;
    int beam_type = 1;
    int ndivx = 200;
    int ndivy = 200;

    int beam_types[] = {1};

    //int attn_idxs[] = {2,3,4,53};
    //int attn_idxs[] = {54};
    int attn_idxs[] = {53};
    //int attn_idxs[] = {3};
    //int attn_idxs[] = {2};

    for (auto beam_type : beam_types)
    {
        vector<double> hole_diameter_array = {0.14,0.16,0.18,0.2,0.22,0.24};
        //vector<double> beam_radius_array = {5,7.5,10,15};
        vector<double> beam_radius_array = {5,10,12.5,15};

        if (beam_type==0) {}
        if (beam_type==1) {
            //ndivx = 20; ndivy = 20;
            //hole_diameter_array.clear();
            //hole_diameter_array.push_back(0.14);
            //hole_diameter_array.push_back(0.16);
            //hole_diameter_array.push_back(0.18);
            //beam_radius_array.clear();
            //beam_radius_array.push_back(8);
        }

        for (auto attn_idx : attn_idxs)
        {
            double wActive = 36; // width of active tantalium aread
            double hActive = 36; // height of active tantalium aread
            double rActive = 0.5*36; // height of active tantalium aread
            int rmdr = 1;
            double dSampleCount = 0.05;
            double attnResolutionMax = 0.03;

            //vector<double> beam_radius_array = {10};
            double attn_factor = 0;
            int attn_constant = 1;
            int attn_exponent = 0;
            if (attn_idx==1) { attn_exponent = 1; attn_factor = 1.E-1; }
            else if (attn_idx==2) { attn_exponent = 2; attn_factor = 1.E-2; }
            else if (attn_idx==3) { attn_exponent = 3; attn_factor = 1.E-3; dSampleCount = 0.2; attnResolutionMax = 0.1;}
            else if (attn_idx==4) { attn_exponent = 4; attn_factor = 1.E-4; dSampleCount = 0.5; attnResolutionMax = 0.8;
                rmdr = 0;
                //hole_diameter_array.clear();
                //hole_diameter_array.push_back(0.14); 
                //hole_diameter_array.push_back(0.16); 
                //hole_diameter_array.push_back(0.18); 
            }
            else if (attn_idx==53) { attn_constant = 5; attn_exponent = 3; attn_factor = attn_constant*1.E-3; dSampleCount = 0.2; attnResolutionMax = 0.2; }
            else if (attn_idx==43) { attn_constant = 4; attn_exponent = 3; attn_factor = attn_constant*1.E-3; dSampleCount = 0.2; attnResolutionMax = 0.2; }
            else if (attn_idx==54) { attn_constant = 5; attn_exponent = 4; attn_factor = attn_constant*1.E-4; dSampleCount = 0.2; attnResolutionMax = 0.2; }

            if (beam_type==1) {
                dSampleCount = dSampleCount*0.5;
            }

            auto drawS = new LKDrawing("drawSummary"); 
            auto frame = new TH2D("frame",";sample radius (mm); attenuation resolution(?)",100,0.1,0.3,100,0,attnResolutionMax);
            drawS -> Add(frame);

            auto nSamples = beam_radius_array.size();
            for (auto iSample=0; iSample<nSamples; ++iSample)
            {
                double beam_radius = beam_radius_array[iSample];
                int beam_size_number = int(beam_radius);
                double beam_sigma = beam_radius/3.;

                TString cname = Form("cvs_%dEm%d_%d",attn_constant,attn_exponent,beam_size_number);
                auto nTestHoles = hole_diameter_array.size();
                auto cvs = LKPainter::GetPainter() -> CanvasResize(cname,600*nTestHoles,1800);
                cvs -> Divide(nTestHoles,3,0.001,0.001);

                auto graphSDV = NewGraphErrors(Form("graphS_%dEm%d_%d",attn_constant,attn_exponent,beam_size_number));
                graphSDV -> SetLineColor(iSample+1);
                drawS -> Add(graphSDV);

                for (auto iHole=0; iHole<nTestHoles; ++iHole)
                {
                    double hole_diameter = hole_diameter_array[iHole];
                    //double hole_diameter = 0.14; // diameter of the hole
                    double hole_radius = 0.5*hole_diameter; // radius of the hole

                    double x0 = 0.5*(wFull - wActive);
                    double y0 = 0.5*(wFull - hActive);
                    double hole_size = hole_radius*hole_radius*TMath::Pi();
                    double attn_size = wActive*hActive;
                    int hole_number = int(100*hole_diameter);

                    double number_of_holes = (attn_size*attn_factor)/hole_size;
                    double x_dist = sqrt(2*hole_size/sqrt(3)/attn_factor);
                    double y_dist = sqrt(3)*x_dist/2;

                    auto ny = int(hActive/y_dist)+1;
                    auto nx1 = int(wActive/x_dist)+1;
                    auto nx = nx1;

                    //cvs -> cd(3*iHole+1);
                    auto draw1 = new LKDrawing();
                    draw1 -> SetCanvas(cvs -> cd(nTestHoles*0+iHole+1));

                    auto draw2 = new LKDrawing("draw2");
                    draw2 -> SetCanvas(cvs -> cd(nTestHoles*1+iHole+1));
                    draw2 -> SetCanvasMargin(0.115,0.14,0.125,0.13);

                    auto draw3 = new LKDrawing("draw3");
                    draw3 -> SetCanvas(cvs -> cd(nTestHoles*2+iHole+1));
                    draw3 -> SetCanvasMargin(0.12,0.12,0.12,0.12);

                    TString hname = Form("hist_%dEm%d_h%d_s%d",attn_constant,attn_exponent,hole_number,beam_size_number);
                    auto hist = new TH2D(hname,";x (mm);y (mm)",100,0,wFull,100,0,wFull);
                    hist -> SetStats(0);
                    draw1 -> Add(hist);

                    vector<TGraph*> graphHoleArray;
                    TGraph* graphHolesInOne;
                    bool useSimpleHoles = (hole_radius<1);
                    if (useSimpleHoles) {
                        graphHolesInOne = NewGraph("graphHolesInOne", 24, (attn_exponent<3?0.4:0.6));
                        graphHoleArray.push_back(graphHolesInOne);
                    }

                    auto graphActive = NewGraph("graphActive");
                    graphActive -> SetLineStyle(2);
                    for (auto i=0; i<=100; ++i)
                        graphActive -> SetPoint(graphActive->GetN(), rActive*cos(i*TMath::Pi()/50)+xMid, rActive*sin(i*TMath::Pi()/50)+xMid);
                    draw1 -> Add(graphActive,"samel");

                    for (auto yy : {y0,xMid,wFull-y0}) {
                        auto line = new TLine(0,yy,wFull,yy);
                        line -> SetLineStyle(2);
                        draw1 -> Add(line,"samel");
                    }
                    for (auto xx : {x0,xMid,wFull-x0}) {
                        auto line = new TLine(xx,0,xx,wFull);
                        line -> SetLineStyle(2);
                        draw1 -> Add(line,"samel");
                    }

                    double xLow  = x0;
                    double xHigh = x0 + (nx-1)*x_dist;
                    double yLow  = y0;
                    double yHigh = y0 + (ny-1)*y_dist;
                    double dxi = (xMid-0.5*(xLow+xHigh));
                    double dyi = (xMid-0.5*(yLow+yHigh));
                    double xi = x0+dxi;
                    double yi = y0+dyi;
                    vector<TVector3> points;
                    for (auto iy=0; iy<ny; ++iy)
                    {
                        nx = nx1;
                        //if (iy%2==1) nx = nx2;
                        double yh = yi + iy*y_dist;
                        for (auto ix=0; ix<nx; ++ix)
                        {
                            double xh = xi + ix*x_dist;
                            if (iy%2==rmdr) xh = xh + 0.5*x_dist;
                            if (xh>x0+wActive) continue;
                            if ((xMid-xh)*(xMid-xh)+(xMid-yh)*(xMid-yh)>rActive*rActive) continue;
                            if (useSimpleHoles) {
                                graphHolesInOne -> SetPoint(graphHolesInOne->GetN(), xh, yh);
                            }
                            else {
                                auto graphHole = NewGraph("graphHoles_%d_%d",ix,iy);
                                graphHole -> SetLineColor(kCyan+2);
                                graphHoleArray.push_back(graphHole);
                                for (auto i=0; i<=20; ++i)
                                    graphHole -> SetPoint(graphHole->GetN(), hole_radius*cos(i*TMath::Pi()/10)+xh, hole_radius*sin(i*TMath::Pi()/10)+yh);
                            }
                            points.push_back(TVector3(xh,yh,0));
                        }
                    }
                    for (auto graphHole : graphHoleArray)
                        draw1 -> Add(graphHole,(useSimpleHoles?"samep":"samel"));

                    TString title = Form("%dx10^{-%d}, #phi_{H}=%.2f, n_{H}=%d, dist_{H}=%.2f, r_{S}=%d",attn_constant,attn_exponent,hole_diameter,int(points.size()),x_dist,beam_size_number);
                    hist -> SetTitle(title);

                    if (1)
                    {
                        auto graphCircle0 = NewGraph("graphCircle0",20,0.6,kBlue);
                        auto graphCircle3 = NewGraph("graphCircle3"); graphCircle3 -> SetLineColor(kBlue);
                        auto graphCircle2 = NewGraph("graphCircle2"); graphCircle2 -> SetLineColor(kBlue);
                        auto graphCircle1 = NewGraph("graphCircle1"); graphCircle1 -> SetLineColor(kBlue);
                        auto graphSimCenter = NewGraph("graphSimCenter");
                        auto graphSimCenterRange = NewGraph("graphSimCenterRange");
                        //double xs1 = xMid - 0.5*x_dist;
                        //double xs2 = xMid + 0.5*x_dist;
                        //double ys1 = xMid - 0.5*y_dist;
                        //double ys2 = xMid + 0.5*y_dist;
                        double xs1 = xMid - 1.0*x_dist;
                        double xs2 = xMid + 1.0*x_dist;
                        double ys1 = xMid - 1.0*y_dist;
                        double ys2 = xMid + 1.0*y_dist;
                        graphSimCenterRange -> SetPoint(0,xs1,ys1);
                        graphSimCenterRange -> SetPoint(1,xs1,ys2);
                        graphSimCenterRange -> SetPoint(2,xs2,ys2);
                        graphSimCenterRange -> SetPoint(3,xs2,ys1);
                        graphSimCenterRange -> SetPoint(4,xs1,ys1);
                        int count = 0;
                        auto histSampleCount1 = new TH1D(hname+"_1d",title+Form(" [Attn./%dx10^{-%d}]",attn_constant,attn_exponent)+Form(";Attenuation / %dx10^{-%d}",attn_constant,attn_exponent),100,1-dSampleCount,1+dSampleCount);
                        auto histSampleCount2 = new TH2D(hname+"_2d",title+Form(" [Attn./%dx10^{-%d}]",attn_constant,attn_exponent)+";sample-center-x (mm); sample-center-y (mm)",ndivx,xs1,xs2,ndivy,ys1,ys2);
                        histSampleCount2 -> SetContour(200);
                        histSampleCount2 -> SetStats(0);
                        histSampleCount2 -> GetZaxis() -> SetMaxDigits(6);
                        histSampleCount2 -> GetZaxis() -> SetDecimals(1);
                        for (auto ix=0; ix<ndivx; ++ix)
                        {
                            //if (beam_type==1) cout << "b=" << beam_type << " " << ", a=" << attn_factor << ", r_b=" << beam_radius << ", d_h=" << hole_diameter << ", x = " << ix << " / " << ndivx << endl;
                            double x1 = xs1 + ix*(xs2-xs1)/ndivx;
                            for (auto iy=0; iy<ndivy; ++iy)
                            {
                                double y1 = ys1 + iy*(ys2-ys1)/ndivy;

                                double total = 0;
                                for (auto point : points) {
                                    double x2 = point.x();
                                    double y2 = point.y();
                                    if ( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) > beam_radius*beam_radius)
                                        continue;
                                    double value = 0;
                                    if (beam_type==0) value = CircleIntersectionArea(x1, y1, beam_radius, x2, y2, hole_radius);
                                    if (beam_type==1) value = ApproximateGaussianIntegral(x1, y1, beam_sigma, beam_sigma, x2, y2, hole_radius);
                                    total += value;
                                }
                                if (beam_type==0) total = total / (TMath::Pi()*beam_radius*beam_radius);
                                if (beam_type==1) total = total / 1;
                                double countRatio = total/attn_factor;
                                histSampleCount2 -> SetBinContent(ix+1,iy+1,countRatio);
                                histSampleCount1 -> Fill(countRatio);
                                //cout << setw(3) << count << setw(10) << x1 << setw(10) << y1 << setw(15) << total << endl;

                                if (ix==0 && iy==0)
                                {
                                    graphCircle0 -> SetPoint(graphCircle0->GetN(),x1,y1);
                                    for (auto i=0; i<=100; ++i) {
                                        graphCircle3 -> SetPoint(graphCircle3->GetN(), beam_radius*3/3*cos(i*TMath::Pi()/50)+x1, beam_radius*3/3*sin(i*TMath::Pi()/50)+y1);
                                        graphCircle2 -> SetPoint(graphCircle2->GetN(), beam_radius*2/3*cos(i*TMath::Pi()/50)+x1, beam_radius*2/3*sin(i*TMath::Pi()/50)+y1);
                                        graphCircle1 -> SetPoint(graphCircle1->GetN(), beam_radius*1/3*cos(i*TMath::Pi()/50)+x1, beam_radius*1/3*sin(i*TMath::Pi()/50)+y1);
                                    }
                                }
                                else
                                    graphSimCenter -> SetPoint(graphSimCenter->GetN(),x1,y1);
                                count++;
                            }
                        }
                        graphSimCenterRange -> SetLineColor(kGreen);
                        draw1 -> Add(graphSimCenterRange,"samel");
                        graphSimCenter -> SetMarkerStyle(20);
                        graphSimCenter -> SetMarkerSize(0.6);
                        graphSimCenter -> SetMarkerColor(kRed);
                        //graphSimCenterRange -> Draw("samel");
                        //graphCircle0 -> Draw("samel");
                        //graphCircle3 -> Draw("samel");
                        //graphCircle0 -> Draw("samep");
                        //graphCircle3 -> Draw("samep");
                        //graphSimCenter -> Draw("samep");
                        draw1 -> Add(graphCircle0,"samel");
                        draw1 -> Add(graphCircle3,"samel");
                        if (beam_type==1) draw1 -> Add(graphCircle2,"samel");
                        if (beam_type==1) draw1 -> Add(graphCircle1,"samel");

                        cvs -> cd(nTestHoles*1+iHole+1);
                        //cvs -> cd(3*iHole+2);
                        draw2 -> Add(histSampleCount2,"colz");
                        auto sdv = histSampleCount1 -> GetStdDev();
                        graphSDV -> SetPoint(graphSDV->GetN(),hole_diameter,sdv);

                        //cvs -> cd(3*iHole+3);
                        draw3 -> Add(histSampleCount1);
                        //histSampleCount1 -> Draw();

                        //cvsCount -> cd(iHole+1);
                        //histSampleCount1 -> DrawClone();
                    }

                    draw1 -> Draw();
                    draw2 -> Draw();
                    draw3 -> SetPaveDx(0.8);
                    draw3 -> SetPaveLineDy(0.1);
                    draw3 -> SetStatsFillStyle(0);
                    draw3 -> Draw();
                }
                cvs -> SaveAs(Form("figures/%s.png",cvs->GetName()));
            }
            //TString cnameS = Form("cvsSummary_%dEm%d",attn_constant,attn_exponent);
            //auto cvsS = LKPainter::GetPainter() -> Canvas(cnameS);
            //drawS -> SetCanvas(cvsS);
            //drawS -> Draw();
            //cvsS -> SaveAs(Form("figures/%s.png",cvsS->GetName()));
        }
    }
}

double CircleIntersectionArea(double x1, double y1, double r1, double x2, double y2, double r2)
{
    double d = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    if (d >= r1 + r2) {
        return 0.0;
    }
    if (d <= std::abs(r1 - r2)) {
        double smallerRadius = std::min(r1, r2);
        return M_PI * smallerRadius * smallerRadius;
    }
    double r1Squared = r1 * r1;
    double r2Squared = r2 * r2;
    double angle1 = std::acos((d * d + r1Squared - r2Squared) / (2 * d * r1));
    double angle2 = std::acos((d * d + r2Squared - r1Squared) / (2 * d * r2));
    double triangleArea = 0.5 * std::sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2));
    double segmentArea1 = r1Squared * angle1;
    double segmentArea2 = r2Squared * angle2;
    return segmentArea1 + segmentArea2 - triangleArea;
}

double Gaussian2D(double x, double y, double xc, double yc, double sigmaX, double sigmaY) {
    double norm = 1.0 / (2 * TMath::Pi() * sigmaX * sigmaY);
    double dx = (x - xc) / sigmaX;
    double dy = (y - yc) / sigmaY;
    return norm * TMath::Exp(-0.5 * (dx * dx + dy * dy));
}

//double IntegrateGaussian2D(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints)
//{
//    double integral = 0.0;
//    double dTheta = 2 * TMath::Pi() / nPoints; // Step size in angle
//    double dr = r0 / nPoints; // Step size in radius
//
//    for (int i = 0; i < nPoints; ++i) {
//        double r = i * dr;
//        for (int j = 0; j < nPoints; ++j) {
//            double theta = j * dTheta;
//            double x = x0 + r * TMath::Cos(theta);
//            double y = y0 + r * TMath::Sin(theta);
//
//            // Accumulate the Gaussian value at this point
//            integral += Gaussian2D(x, y, xc, yc, sigmaX, sigmaY) * r * dr * dTheta;
//        }
//    }
//
//    return integral;
//}

double IntegrateGaussian2D(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints)
{
    double cx = TMath::Sqrt((x0-xc)*(x0-xc) + (y0-yc)*(y0-yc));
    return IntegrateGaussian2D_SumGamma(sigmaX, cx, r0);
}

double IntegrateGaussian2D_SumGamma(double sigma, double cx, double r0)
{
    lk_debug << "cx = " << cx << endl;
    cx = cx/sigma;
    r0 = r0/sigma;
    double valueSum = 0;
    for (int k=0; k<4; ++k)
    {
        double k_factorial = TMath::Factorial(k);
        double value_add = TMath::Power(cx*cx/2,k)/(k_factorial*k_factorial) * TMath::Gamma(k+1,r0*r0/2);
        lk_debug << k << ") " << TMath::Power(cx*cx/2,k) << " / " << (k_factorial*k_factorial) << " * " << TMath::Gamma(k+1,r0*r0/2) << " = " << value_add << endl;
        valueSum += value_add;
    }
    lk_debug << TMath::Exp(-cx*cx/2) << endl;
    lk_debug << TMath::Exp(-cx*cx/2) * valueSum << endl;
    double integral = 1 - TMath::Exp(-cx*cx/2) * valueSum;
    lk_debug << integral << endl;
    return integral;
}

double ApproximateGaussianIntegral(double xc, double yc, double sigmaX, double sigmaY,
                                   double x0, double y0, double r0) {
    double dx = (x0 - xc) / sigmaX;
    double dy = (y0 - yc) / sigmaY;
    double gaussianValue = (1.0 / (2 * TMath::Pi() * sigmaX * sigmaY)) *
                           TMath::Exp(-0.5 * (dx * dx + dy * dy));

    double circleArea = TMath::Pi() * r0 * r0;

    return circleArea * gaussianValue;
}

TGraph *NewGraph(TString name, int mst, double msz, int mcl, int lst, int lsz, int lcl)
{
    auto graph = new TGraph();
    graph -> SetName(name);
    if (mst<=0) mst = 20;
    if (msz<=0) msz = 1;
    if (mcl<0) mcl = kBlack;
    graph -> SetMarkerStyle(mst);
    graph -> SetMarkerSize(msz);
    graph -> SetMarkerColor(mcl);
    if (lst<0) lst = 1;
    if (lsz<0) lsz = 1;
    if (lcl<0) lcl = mcl;
    graph -> SetLineStyle(lst);
    graph -> SetLineWidth(lsz);
    graph -> SetLineColor(lcl);
    graph -> SetFillStyle(0);
    return graph;
}

TGraphErrors *NewGraphErrors(TString name, int mst, double msz, int mcl, int lst, int lsz, int lcl)
{
    auto graph = new TGraphErrors();
    graph -> SetName(name);
    if (mst<=0) mst = 20;
    if (msz<=0) msz = 1;
    if (mcl<0) mcl = kBlack;
    graph -> SetMarkerStyle(mst);
    graph -> SetMarkerSize(msz);
    graph -> SetMarkerColor(mcl);
    if (lst<0) lst = 1;
    if (lsz<0) lsz = 1;
    if (lcl<0) lcl = mcl;
    graph -> SetLineStyle(lst);
    graph -> SetLineWidth(lsz);
    graph -> SetLineColor(lcl);
    graph -> SetFillStyle(0);
    return graph;
}
