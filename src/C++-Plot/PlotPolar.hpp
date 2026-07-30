#pragma once
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstddef>
#include <format>

template<typename T>
using Vec  = std::vector<T>;
using dVec = Vec<double>;
using Point = std::pair<double, double>;

constexpr double PI = 3.141592653589793; // \pi

double estimateStringWidth(const std::string& text, double fontsize)
{
    std::string text_ = text;
    size_t pos = 0;

    // --- First loop: Replace HTML entities like &#XXXX; with '^' ---
    while (true)
    {
        size_t startPos = text_.find("&#", pos);
        if (startPos == std::string::npos) {
            break; // No more "&#" sequences found, exit loop
        }

        size_t endPos = text_.find(";", startPos); // Search for ';' AFTER the found "&#"
        if (endPos != std::string::npos)
        {
            // Found a complete "&#...;" sequence
            text_.replace(startPos, endPos - startPos + 1, "^");
            // Advance 'pos' to continue searching immediately after the replacement
            pos = startPos + 1; // '^' is 1 character long
        }
        else
        {
            // Found "&#" but no matching ";". To prevent an infinite loop,
            // advance 'pos' past the found "&#" and continue searching.
            pos = startPos + 2;
        }
    }

    // --- Second loop: Replace "&amp;" with '&' ---
    pos = 0; // Reset position for the new search
    while (true)
    {
        size_t startPos = text_.find("&amp;", pos);
        if (startPos == std::string::npos) {
            break; // No more "&amp;" sequences found, exit loop
        }
        // Found "&amp;"
        text_.replace(startPos, 5, "&"); // " &amp;" is 5 characters long
        // Advance 'pos' to continue searching immediately after the replacement
        pos = startPos + 1; // '&' is 1 character long
    }

    // --- Character width estimation (this part seems fine) ---
    double totalRelativeWidth = 0.0;
    for ( auto c : text_ )
    {
        if (std::strchr("lij|' !.,;:()[]{}",c)!=nullptr) totalRelativeWidth += 0.3;
        else if (std::strchr("fzrstc-+^",c)!=nullptr) totalRelativeWidth += 0.4;
        else if (std::strchr("WwMmOQG%&",c)!=nullptr) totalRelativeWidth += 0.9;
        else if (std::isupper(static_cast<unsigned char>(c))) totalRelativeWidth += 0.7;
        else totalRelativeWidth += 0.6;
    }
    return totalRelativeWidth*fontsize;
}

enum class LegendPos
{
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight
};

class PlotPolarData
{
    public:
        dVec rs;
        dVec thetas;
        bool pointsFill=false;
        double pointsBorderWidth=1.0;
        double pointsRadius=5.0;
        double pointsBorderOpacity=1.0;
        double pointsFillOpacity=1.0;
        std::string pointsBorderColor="#000000";
        std::string pointsFillColor="#0000ff";
        std::string pointsLabel="";
        PlotPolarData(const dVec& r,const dVec& t, bool pf=false, double pbw=1.0, double pr=5.0, double pbo=1.0, double pfo=1.0, size_t pbcr=0, size_t pbcg=0, size_t pbcb=0, size_t pfcr=0, size_t pfcg=0, size_t pfcb=255,  std::string label="")
        : rs(r), thetas(t), pointsFill(pf), pointsBorderWidth(pbw), pointsRadius(pr), pointsBorderOpacity(pbo), pointsFillOpacity(pfo), pointsLabel(label)
        {
            // Border Color
            if ( (pbcr<256) && (pbcg<256) && (pbcb<256) )
            {
                pointsBorderColor = "#"+((pbcr<16) ? "0"+std::format("{:x}",pbcr) : std::format("{:x}",pbcr))
                    +((pbcg<16) ? "0"+std::format("{:x}",pbcg) : std::format("{:x}",pbcg))
                    +((pbcb<16) ? "0"+std::format("{:x}",pbcb) : std::format("{:x}",pbcb));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to black points' border...";
            }
            // Fill Color
            if ( (pbcr<256) && (pbcg<256) && (pbcb<256) )
            {
                pointsFillColor = "#"+((pfcr<16) ? "0"+std::format("{:x}",pfcr) : std::format("{:x}",pfcr))
                    +((pfcg<16) ? "0"+std::format("{:x}",pfcg) : std::format("{:x}",pfcg))
                    +((pfcb<16) ? "0"+std::format("{:x}",pfcb) : std::format("{:x}",pfcb));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to blue points' fill...";
            }
        }
};

class PlotPolar
{
    private:
        Vec<PlotPolarData> plotData;
        double epsilon = 1.0e-10;
        double width = 1000.0;
        dVec rVals;
        dVec thetaVals;
        double padding = 100.0;
        double plotPad = 10.0;
        double rMax;
        double rMin;
        double rRange;
        double rScale;
        dVec xVals;
        dVec yVals;
        double drawW;
        size_t numRDivisions = 5.0;
        size_t numThetaDivisions = 12;
        dVec divisionsR;
        dVec divisionsTheta;
        double rStep;
        double thetaStep;
        double gLeft;
        double gRight;
        double gTop;
        double gBottom;
        bool border=true;
        double borderLineWidth = 2.5;
        std::string borderLineColor = "#000000";
        // double axisLineWidth = 1.5;
        // std::string axisLineColor = "#000000";
        double gridLineWidth = 0.5;
        std::string gridLineColor = "#101010";
        bool plotNumerals = false;
        size_t plotNumeralsSignificantDigits = 3;
        double plotNumeralsFontSize = 15;
        std::string plotNumeralsColor = "#3f003f";
        bool plotTitle = false;
        double plotTitleFontSize = 36;
        std::string plotTitleText = "Title!";
        std::string plotTitleColor = "#000000";
        bool plotPointsFill = true;
        double plotPointsBorderWidth = 1.5;
        double plotPointsRadius = 5;
        double plotPointsBorderOpacity = 1.0;
        double plotPointsFillOpacity = 1.0;
        std::string plotPointsBorderColor = "#000000";
        std::string plotPointsFillColor = "#0000ff";
        std::string plotPointsLabel = "";
        bool plotLegend = false;
        LegendPos legendPos = LegendPos::TopLeft;
        double legendFontSize = 12;
        std::string backColor = "#ffffff";
        inline double calcRDivisionsStep()
        {
            double rawStep  = rRange / numRDivisions;
            double exponent = std::floor(std::log10(rawStep));
            double fraction = rawStep / std::pow(10.0,exponent);
            double cleanFraction = 1.0;
            if (fraction<1.2) cleanFraction = 1.0;
            else if (fraction<1.5) cleanFraction = 1.25;
            else if (fraction<2.25) cleanFraction = 2.0;
            else if (fraction<3.75) cleanFraction = 2.5;
            else if (fraction<7.5) cleanFraction = 5.0;
            else cleanFraction = 10.0;
            return cleanFraction*std::pow(10.0,exponent);
        };
        inline double toPixelX(double x) { return (width/2.0 + x*rScale); };
        inline double toPixelY(double y) { return (width/2.0 - y*rScale); };
        inline double rttox(double r, double t) { return r*cos(t); };
        inline double rttoy(double r, double t) { return r*sin(t); };
        inline void writeNumerals(std::ofstream& file);
        inline void writeplotTitle(std::ofstream& file);
        inline void drawLegend(std::ofstream& file);
    public:
        PlotPolar(const Vec<PlotPolarData>& data, double w=1000.0, double pad=100.0, double ppad = 50.0, size_t nrd = 5, size_t ntd = 12)
        : width(w), padding(pad), plotPad(ppad), numRDivisions(nrd), numThetaDivisions(ntd)
        {
            if ( width<=2*padding )
                throw std::runtime_error("Padding is too large! It exceeds width!");
            for ( auto pd : data )
            {
                rRange = 0.0;
                if ( rVals.empty() )
                    throw std::runtime_error("Values of Rs are empty!");
                if ( rVals.size() != thetaVals.size() )
                    throw std::runtime_error("Rs, and Thetas don't match in size!");
                double rMinLocal = *std::min_element(pd.rs.begin(),pd.rs.end());
                double rMaxLocal = *std::max_element(pd.rs.begin(),pd.rs.end());
                double rRangeLocal = (std::abs(rMinLocal)<=std::abs(rMaxLocal)) ? std::abs(rMaxLocal) : std::abs(rMinLocal);
                rRange = (rRangeLocal<=rRange) ? rRange : rRangeLocal;
                plotData.push_back(pd);
            }
            drawW = width - 2*padding;
            rScale = 0.5 * drawW / rRange;
        }
        PlotPolar(dVec rs, dVec thetas, double w = 1000.0, double pad = 100.0, double ppad = 10.0, size_t nrd = 5, size_t ntd = 12)
        : rVals(rs), thetaVals(thetas), width(w), padding(pad), plotPad(ppad), numRDivisions(nrd), numThetaDivisions(ntd)
        {
            if ( rVals.empty() )
                throw std::runtime_error("Values of Rs are empty!");
            if ( rVals.size() != thetaVals.size() )
                throw std::runtime_error("Rs, and Thetas don't match in size!");
            if ( width<=2*padding )
                throw std::runtime_error("Padding is too large! It exceeds width!");
            rMin = *std::min_element(rVals.begin(),rVals.end());
            rMax = *std::max_element(rVals.begin(),rVals.end());
            rRange = (std::abs(rMin)<=std::abs(rMax)) ? std::abs(rMax) : std::abs(rMin);
            drawW = width - 2*padding;
            rScale = 0.5 * drawW / rRange;
        };
        inline void addData(const dVec& r, const dVec& t, bool pf, double pbw=1.0, double pr=5.0, double pbo=1.0, double pfo=1.0, double pbcr=0, double pbcg=0, double pbcb=0, double pfcr=0, double pfcg=0, double pfcb=255, std::string label="")
        {
            if ( r.empty() )
                throw std::runtime_error("Values of Rs are empty!");
            if ( r.size() != t.size() )
                throw std::runtime_error("Rs, and Thetas don't match in size!");
            PlotPolarData ppd(r,t,pf,pbw,pr,pbo,pfo,pbcr,pbcg,pbcb,pfcr,pfcg,pfcb,label);
            plotData.push_back(ppd);
            double rMinLocal = *std::min_element(r.begin(),r.end());
            double rMaxLocal = *std::max_element(r.begin(),r.end());
            double rRangeLocal = (std::abs(rMinLocal)<=std::abs(rMaxLocal)) ? std::abs(rMaxLocal) : std::abs(rMinLocal);
            rRange = (rRangeLocal<=rRange) ? rRange : rRangeLocal;
            rScale = 0.5 * drawW / rRange;
        }
        inline void addMultipleData(const Vec<PlotPolarData>& data)
        {
            for ( auto pd : data )
            {
                if ( pd.rs.empty() )
                    throw std::runtime_error("Values of Rs are empty!");
                if ( pd.rs.size() != pd.thetas.size() )
                    throw std::runtime_error("Rs, and Thetas don't match in size!");
                plotData.push_back(pd);
                double rMinLocal = *std::min_element(pd.rs.begin(),pd.rs.end());
                double rMaxLocal = *std::max_element(pd.rs.begin(),pd.rs.end());
                double rRangeLocal = (std::abs(rMinLocal)<=std::abs(rMaxLocal)) ? std::abs(rMaxLocal) : std::abs(rMinLocal);
                rRange = (rRangeLocal<=rRange) ? rRange : rRangeLocal;
            }
            rScale = 0.5 * drawW / rRange;
        }
        inline void forceRMax(double rm)
        {
            rMax = std::abs(rm);
            rRange = rMax;
            rScale = 0.5 * drawW / rRange;
        }

        inline void plotSVG(std::string filename);

        inline void setnumRDivisions(size_t nrd=5) {numRDivisions=nrd;};
        inline void setnumThetaDivisions(size_t ntd=8) {numThetaDivisions=ntd;};
        inline void setplotPad(double pp) {plotPad=pp;};
        inline void setbackColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                backColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to white background...";
            }
        };
        inline void setplotTitle(bool title_) {plotTitle=title_;};
        inline void setplotTitleFontSize(double ptfs) {plotTitleFontSize=ptfs;};
        inline void setplotTitleText(std::string titletext) {plotTitleText=titletext;};
        inline void setplotTitleColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                plotTitleColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to black plot title...";
            }
        };
        inline void setborder(bool border_) {border=border_;};
        inline void setborderLineWidth(double blw) {borderLineWidth=blw;};
        inline void setborderLineColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                borderLineColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to black borders...";
            }
        };
        // inline void setaxisLineWidth(double alw) {axisLineWidth=alw;};
        // inline void setaxisLineColor(size_t R, size_t G, size_t B)
        // {
        //     if ( (R<256) && (G<256) && (B<256) )
        //     {
        //         axisLineColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
        //             +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
        //             +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
        //     }
        //     else
        //     {
        //         std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to black axes lines...";
        //     }
        // };
        inline void setgridLineWidth(double glw) {gridLineWidth=glw;};
        inline void setgridLineColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                gridLineColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to gray grid lines...";
            }
        };
        inline void setplotNumerals(bool pn) {plotNumerals=pn;};
        inline void setplotNumeralsFontSize(double pnfs) {plotNumeralsFontSize=pnfs;};
        inline void setplotNumeralsSignificantDigits(size_t pnsd) {plotNumeralsSignificantDigits=pnsd;};
        inline void setplotNumeralsColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                plotNumeralsColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to black numerals...";
            }
        };
        inline void setplotPointsFill(bool ppf) {plotPointsFill=ppf;};
        inline void setplotPointsBorderWidth(double ppbw) {plotPointsBorderWidth=ppbw;};
        inline void setplotPointsRadius(double ppr) {plotPointsRadius=ppr;};
        inline void setplotPointsBorderOpacity(double op) {plotPointsBorderOpacity=op;};
        inline void setplotPointsFillOpacity(double op) {plotPointsFillOpacity=op;};
        inline void setplotPointsBorderColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                plotPointsBorderColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to black plot points border...";
            }
        };
        inline void setplotPointsFillColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                plotPointsFillColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to blue plot points fill...";
            }
        };
        inline void setplotPointsLabel(std::string ppl) {plotPointsLabel=ppl;};
        inline void setplotLegend(bool pl) {plotLegend=pl;};
        inline void setlegendPos(LegendPos lp) {legendPos=lp;};
        inline void setlegendFontSize(double lfs) {legendFontSize=lfs;};
};

inline void PlotPolar::writeNumerals(std::ofstream& file)
{
    file << " \n";
    file << " <g font-family=\"sans-serif\" font-size=\"" << PlotPolar::plotNumeralsFontSize
         << "\" fill=\"" << PlotPolar::plotNumeralsColor << "\">\n";

    // --- X-AXIS NUMERALS ---
    for ( size_t i=0; i<PlotPolar::divisionsR.size(); ++i)
    {
        std::stringstream ss;
        // defaultfloat strips trailing zeros and manages scientific boundaries dynamically
        // 3 significant digits means: 1234, 123.4, 12.34, 1.234, or 1.234e+05
        ss << std::defaultfloat << std::setprecision(PlotPolar::plotNumeralsSignificantDigits) << PlotPolar::divisionsR[i];
        std::string labelStr = ss.str();

        // Elegant parsing check for scientific strings
        size_t ePos = labelStr.find('e');
        if (ePos != std::string::npos)
        {
            std::string mantissa = labelStr.substr(0, ePos);
            std::string exponent = labelStr.substr(ePos + 1);
            int expVal = std::stoi(exponent); // Strips leading '+' or zeros (e.g. "+04" -> 4)
            std::string superscriptExp = "";
            std::string expStr = std::to_string(expVal);

            for (char c : expStr)
            {
                if (c == '-') superscriptExp += "&#8315;"; // Unicode ⁻
                else if (c == '0') superscriptExp += "&#8304;";  // Unicode ⁰
                else if (c == '1') superscriptExp += "&#185;";  // Unicode ¹
                else if (c == '2') superscriptExp += "&#178;";  // Unicode ²
                else if (c == '3') superscriptExp += "&#179;";  // Unicode ³
                else if (c == '4') superscriptExp += "&#8308;"; // Unicode ⁴
                else if (c == '5') superscriptExp += "&#8309;"; // Unicode ⁵
                else if (c == '6') superscriptExp += "&#8310;"; // Unicode ⁶
                else if (c == '7') superscriptExp += "&#8311;"; // Unicode ⁷
                else if (c == '8') superscriptExp += "&#8312;"; // Unicode ⁸
                else if (c == '9') superscriptExp += "&#8313;"; // Unicode ⁹
            }
            labelStr = mantissa + "&#215;10" + superscriptExp;
        }

        file << "   <text x=\"" << PlotPolar::toPixelX(PlotPolar::divisionsR[i])
             << "\" y=\"" << PlotPolar::width*0.51 + PlotPolar::plotNumeralsFontSize
             << "\" text-anchor=\"middle\" dominant-baseline=\"hanging\">"
             << labelStr << "</text>\n";
    }
    file << "   <text x=\"" << PlotPolar::toPixelX(0.0)
         << "\" y=\"" << PlotPolar::width*0.51 + PlotPolar::plotNumeralsFontSize
         << "\" text-anchor=\"middle\" dominant-baseline=\"hanging\">"
         << "0" << "</text>\n";

    // --- Y-AXIS NUMERALS ---
    double toDegree = 180/PI;
    for ( size_t i=0; i<PlotPolar::divisionsTheta.size(); ++i )
    {
        std::stringstream ss;
        ss << std::defaultfloat << std::setprecision(PlotPolar::plotNumeralsSignificantDigits) << PlotPolar::divisionsTheta[i]*toDegree;
        std::string labelStr = ss.str();

        // Elegant parsing check for scientific strings
        size_t ePos = labelStr.find('e');
        if (ePos != std::string::npos)
        {
            std::string mantissa = labelStr.substr(0, ePos);
            std::string exponent = labelStr.substr(ePos + 1);
            int expVal = std::stoi(exponent);
            std::string superscriptExp = "";
            std::string expStr = std::to_string(expVal);
            for (char c : expStr)
            {
                if (c == '-') superscriptExp += "&#8315;"; // Unicode ⁻
                else if (c == '0') superscriptExp += "&#8304;";  // Unicode ⁰
                else if (c == '1') superscriptExp += "&#185;";  // Unicode ¹
                else if (c == '2') superscriptExp += "&#178;";  // Unicode ²
                else if (c == '3') superscriptExp += "&#179;";  // Unicode ³
                else if (c == '4') superscriptExp += "&#8308;"; // Unicode ⁴
                else if (c == '5') superscriptExp += "&#8309;"; // Unicode ⁵
                else if (c == '6') superscriptExp += "&#8310;"; // Unicode ⁶
                else if (c == '7') superscriptExp += "&#8311;"; // Unicode ⁷
                else if (c == '8') superscriptExp += "&#8312;"; // Unicode ⁸
                else if (c == '9') superscriptExp += "&#8313;"; // Unicode ⁹
            }
            labelStr = mantissa + "&#215;10" + superscriptExp;
        }

        // Added dominant-baseline="central" to keep numbers aligned with the tick line rows
        file << "   <text x=\"" <<  PlotPolar::toPixelX(PlotPolar::rttox(((drawW+plotPad*0.75)/(drawW))*rRange,divisionsTheta[i]))
             << "\" y=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(((drawW+plotPad*0.75)/(drawW))*rRange,divisionsTheta[i]))+PlotPolar::plotNumeralsFontSize*0.25
             << "\" text-anchor=\"middle\" dominant-baseline=\"central\">"
             << labelStr << "</text>\n";
    }
    file << " </g>\n";
}

inline void PlotPolar::writeplotTitle(std::ofstream& file)
{
    file << " <g font-family=\"sans-serif\" font-size=\"" << PlotPolar::plotTitleFontSize
         << "\" fill=\"" << PlotPolar::plotTitleColor << "\">\n";
    file << "   <text x=\"" << PlotPolar::width / 2.0 << "\" y=\"" << PlotPolar::plotTitleFontSize+10.0
        << "\" text-anchor=\"middle\" dominant-baseline=\"central\" font-weight=\"bold\">"
        << PlotPolar::plotTitleText << "</text>\n";
    file << " </g>\n";
}

inline void PlotPolar::drawLegend(std::ofstream& file)
{
    double legendTopPad = 20.0;
    double legendBottomPad = 20.0;
    double legendLeftPad = 20.0;
    double legendRightPad = 20.0;
    double legendWidth = 0.0;
    double legendHeight = 0.0;
    if (PlotPolar::rVals.empty())
    {
        if (PlotPolar::plotData.empty())
        {
            legendWidth  = 220.0;
            legendHeight = 220.0;
        }
        else
        {
            double currentWidth = 0.0;
            for ( auto pd : PlotPolar::plotData)
            {
                double newWidth = estimateStringWidth(pd.pointsLabel,PlotPolar::legendFontSize);
                if ( newWidth > currentWidth ) currentWidth = newWidth;
            }
            legendWidth  = currentWidth + 40.0;
            legendHeight = plotData.size()*PlotPolar::legendFontSize*2.0 + 20.0;
        }
    }
    else
    {
        if (PlotPolar::plotData.empty())
        {
            legendWidth  = estimateStringWidth(PlotPolar::plotPointsLabel,PlotPolar::legendFontSize) + 40.0;
            legendHeight = PlotPolar::legendFontSize*2.0 + 20.0;
        }
        else
        {
            double currentWidth = estimateStringWidth(PlotPolar::plotPointsLabel,PlotPolar::legendFontSize);
            for ( auto pd : PlotPolar::plotData)
            {
                double newWidth = estimateStringWidth(pd.pointsLabel,PlotPolar::legendFontSize);
                if ( newWidth > currentWidth ) currentWidth = newWidth;
            }
            legendWidth  = currentWidth + 40.0;
            legendHeight = (plotData.size()+1)*PlotPolar::legendFontSize*2.0 + 20.0;
        }
    }

    double legendX = 0.0;
    double legendY = 0.0;
    switch ( PlotPolar::legendPos )
    {
        case LegendPos::TopLeft:
            legendX = PlotPolar::gLeft-PlotPolar::plotPad+legendLeftPad;
            legendY = PlotPolar::gTop-PlotPolar::plotPad+legendTopPad;
            break;
        case LegendPos::TopRight:
            legendX = PlotPolar::gRight+PlotPolar::plotPad-legendRightPad-legendWidth;
            legendY = PlotPolar::gTop-PlotPolar::plotPad+legendTopPad;
            break;
        case LegendPos::BottomLeft:
            legendX = PlotPolar::gLeft-PlotPolar::plotPad+legendLeftPad;
            legendY = PlotPolar::gBottom+PlotPolar::plotPad-legendBottomPad-legendHeight;
            break;
        case LegendPos::BottomRight:
            legendX = PlotPolar::gRight+PlotPolar::plotPad-legendRightPad-legendWidth;
            legendY = PlotPolar::gBottom+PlotPolar::plotPad-legendBottomPad-legendHeight;
            break;
        default:
            legendX = PlotPolar::gLeft-PlotPolar::plotPad+legendLeftPad;
            legendY = PlotPolar::gTop-PlotPolar::plotPad+legendTopPad;
            break;
    }
    double itemY = 0.0;
    // Drawing the empty legend box
    file << "   <rect x=\"" << legendX << "\" y=\"" << legendY
        << "\" width=\"" << legendWidth << "\" height=\"" << legendHeight
        << "\" fill=\"" << PlotPolar::backColor << "\" fill-opacity=\"0.75\" stroke=\"" << PlotPolar::borderLineColor
        << "\" stroke-width=\"1\" rx=\"4\" />\n";
    if ( PlotPolar::rVals.empty() )
    {
        if (PlotPolar::plotData.empty() )
        {
            std::cout << "Nothing to write in legend!";
        }
        else
        {
            size_t counter = 0;
            for ( auto pd : plotData )
            {
                itemY = legendY + counter*PlotPolar::legendFontSize*2.0 + 10.0;
                if ( pd.pointsFill )
                {
                    file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                        << "\" r=\"" << pd.pointsRadius
                        << "\" stroke=\"" << pd.pointsBorderColor
                        << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                        << "\" stroke-width=\"" << pd.pointsBorderWidth
                        << "\" fill=\"" << pd.pointsFillColor
                        << "\" fill-opacity=\"" << pd.pointsFillOpacity
                        << "\" />\n";
                }
                else
                {
                    file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                        << "\" r=\"" << pd.pointsRadius
                        << "\" stroke=\"" << pd.pointsBorderColor
                        << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                        << "\" fill=\"none\" stroke-width=\"" << pd.pointsBorderWidth << "\" />\n";
                }
                file << "   <text x=\"" << legendX+legendWidth/2.0+10.0 << "\" y=\"" << itemY+1.5*PlotPolar::legendFontSize
                    << "\" font-family=\"sans-serif\" font-size=\"" << PlotPolar::legendFontSize
                    << "\" fill=\"" << pd.pointsBorderColor
                    << "\" text-anchor=\"middle\" dominant-baseline=\"central\">"
                    << pd.pointsLabel << "</text>\n";
                counter += 1;
            }
        }
    }
    else
    {
        if (PlotPolar::plotData.empty())
        {
            std::string dashArrayAttr = "";
            itemY = legendY + 10.0;
            if ( PlotPolar::plotPointsFill )
            {
                file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                    << "\" r=\"" << PlotPolar::plotPointsRadius
                    << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                    << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                    << "\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth
                    << "\" fill=\"" << PlotPolar::plotPointsFillColor
                    << "\" fill-opacity=\"" << PlotPolar::plotPointsFillOpacity
                    << "\" />\n";
            }
            else
            {
                file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                    << "\" r=\"" << PlotPolar::plotPointsRadius
                    << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                    << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                    << "\" fill=\"none\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth << "\" />\n";
            }
            file << "   <text x=\"" << legendX+legendWidth/2.0+10.0 << "\" y=\"" << itemY+1.5*PlotPolar::legendFontSize
                << "\" font-family=\"sans-serif\" font-size=\"" << PlotPolar::legendFontSize
                << "\" fill=\"" << PlotPolar::plotPointsBorderColor
                << "\" text-anchor=\"middle\" dominant-baseline=\"central\">"
                << PlotPolar::plotPointsLabel << "</text>\n";
        }
        else
        {
            size_t counter = 0;
            // writing original data
            itemY = legendY + 10.0;
            if ( PlotPolar::plotPointsFill )
            {
                file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                    << "\" r=\"" << PlotPolar::plotPointsRadius
                    << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                    << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                    << "\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth
                    << "\" fill=\"" << PlotPolar::plotPointsFillColor
                    << "\" fill-opacity=\"" << PlotPolar::plotPointsFillOpacity
                    << "\" />\n";
            }
            else
            {
                file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                    << "\" r=\"" << PlotPolar::plotPointsRadius
                    << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                    << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                    << "\" fill=\"none\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth << "\" />\n";
            }
            file << "   <text x=\"" << legendX+legendWidth/2.0+10.0 << "\" y=\"" << itemY+1.5*PlotPolar::legendFontSize
                << "\" font-family=\"sans-serif\" font-size=\"" << PlotPolar::legendFontSize
                << "\" fill=\"" << PlotPolar::plotPointsBorderColor
                << "\" text-anchor=\"middle\" dominant-baseline=\"central\">"
                << PlotPolar::plotPointsLabel << "</text>\n";
            counter += 1;
            // writing added data
            for ( auto pd : plotData )
            {
                itemY = legendY + counter*PlotPolar::legendFontSize*2.0 + 10.0;
                if ( pd.pointsFill )
                {
                    file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                        << "\" r=\"" << pd.pointsRadius
                        << "\" stroke=\"" << pd.pointsBorderColor
                        << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                        << "\" stroke-width=\"" << pd.pointsBorderWidth
                        << "\" fill=\"" << pd.pointsFillColor
                        << "\" fill-opacity=\"" << pd.pointsFillOpacity
                        << "\" />\n";
                }
                else
                {
                    file << " <circle cx=\"" << legendX + 15.0 << "\" cy=\"" << itemY+PlotPolar::legendFontSize
                        << "\" r=\"" << pd.pointsRadius
                        << "\" stroke=\"" << pd.pointsBorderColor
                        << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                        << "\" fill=\"none\" stroke-width=\"" << pd.pointsBorderWidth << "\" />\n";
                }
                file << "   <text x=\"" << legendX+legendWidth/2.0+10.0 << "\" y=\"" << itemY+1.5*PlotPolar::legendFontSize
                    << "\" font-family=\"sans-serif\" font-size=\"" << PlotPolar::legendFontSize
                    << "\" fill=\"" << pd.pointsBorderColor
                    << "\" text-anchor=\"middle\" dominant-baseline=\"central\">"
                    << pd.pointsLabel << "</text>\n";
                counter += 1;
            }
        }
    }
}

inline void PlotPolar::plotSVG(std::string filename)
{
    std::ofstream file(filename, std::ios::trunc);

    // Setting up some parameters
    rStep = calcRDivisionsStep();
    thetaStep = 2*PI/numThetaDivisions;
    for ( double val=rStep; val<=rRange+epsilon; val+=rStep )
    {
        divisionsR.push_back(val);
    }
    for ( size_t i=0; i<numThetaDivisions; ++i )
    {
        divisionsTheta.push_back(i*thetaStep);
    }

    // Background
    file << "<svg width=\"" << PlotPolar::width << "\" height=\"" << PlotPolar::width
        << "\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:" << PlotPolar::backColor <<";\">\n";

    file << " \n";
    file << " \n";
    // Painting background
    file << " <rect width=\"" << PlotPolar::width << "\" height=\"" << PlotPolar::width
        << "\" fill=\"" << PlotPolar::backColor << "\" />\n";
    // Marking corners
    PlotPolar::gLeft   = PlotPolar::toPixelX(PlotPolar::rttox(rRange, PI));
    PlotPolar::gRight  = PlotPolar::toPixelX(PlotPolar::rttox(rRange, 0.0));
    PlotPolar::gTop    = PlotPolar::toPixelY(PlotPolar::rttoy(rRange, PI/2.0));
    PlotPolar::gBottom = PlotPolar::toPixelY(PlotPolar::rttoy(rRange, -PI/2.0));
    // Border
    if ( border )
    {
        file << " \n";
        file << " <rect x=\"" << PlotPolar::gLeft-PlotPolar::plotPad << "\" y=\"" << PlotPolar::gTop-PlotPolar::plotPad
             << "\" width=\"" << (PlotPolar::gRight-PlotPolar::gLeft+2*PlotPolar::plotPad)
             << "\" height=\"" << (PlotPolar::gBottom-PlotPolar::gTop+2*PlotPolar::plotPad)
             << "\" fill=\"none\" stroke=\"" << PlotPolar::borderLineColor << "\" stroke-width=\"" << PlotPolar::borderLineWidth << "\" />\n";
    }
    // Drawing xAxis
    file << " \n";
    double xAxisX1 = PlotPolar::gLeft;
    double xAxisX2 = PlotPolar::gRight;
    double xAxisY1 = PlotPolar::toPixelY(0.0);
    double xAxisY2 = PlotPolar::toPixelY(0.0);
    // file << " <line x1=\"" << xAxisX1 << "\" y1=\"" << xAxisY1
    //     << "\" x2=\"" << xAxisX2 << "\" y2=\"" << xAxisY2
    //     << "\" stroke=\"" << PlotPolar::axisLineColor << "\" stroke-width=\"" << PlotPolar::axisLineWidth << "\" stroke-linecap=\"round\" />\n";
    // Drawing yAxis
    double yAxisX1 = PlotPolar::toPixelX(0.0);
    double yAxisX2 = PlotPolar::toPixelX(0.0);
    double yAxisY1 = PlotPolar::gBottom;
    double yAxisY2 = PlotPolar::gTop;
    // file << " <line x1=\"" << yAxisX1 << "\" y1=\"" << yAxisY1
    //     << "\" x2=\"" << yAxisX2 << "\" y2=\"" << yAxisY2
    //     << "\" stroke=\"" << PlotPolar::axisLineColor << "\" stroke-width=\"" << PlotPolar::axisLineWidth << "\" stroke-linecap=\"round\" />\n";
    // Drawing the grid
    for ( double r : PlotPolar::divisionsR )
    {
        file << " <circle cx=\"" << yAxisX1 << "\" cy=\"" << xAxisY1
            << "\" r=\"" << r*rScale << "\" stroke=\"" << PlotPolar::gridLineColor
            << "\" fill=\"none\" stroke-width=\"" << PlotPolar::gridLineWidth << "\" />\n";
    }
    // file << " <circle cx=\"" << yAxisX1 << "\" cy=\"" << xAxisY1
    //     << "\" r=\"" << divisionsR[divisionsR.size()-1]*rScale << "\" stroke=\"" << PlotPolar::axisLineColor
    //     << "\" fill=\"none\" stroke-width=\"" << PlotPolar::axisLineWidth << "\" />\n";
    for ( double t : divisionsTheta )
    {
        file << " <line x1=\"" << yAxisX1 << "\" y1=\"" << xAxisY1
            << "\" x2=\"" << PlotPolar::toPixelX(PlotPolar::rttox(rRange,t))
            << "\" y2=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(rRange,t))
            << "\" stroke=\"" << PlotPolar::gridLineColor
            << "\" stroke-width=\"" << PlotPolar::gridLineWidth << "\" stroke-linecap=\"round\" />\n";
    }
    // Writing Numbers
    if ( PlotPolar::plotNumerals )
    {
        PlotPolar::writeNumerals(file);
    }

    // Plot
    if ( rVals.empty() )
    {
        if ( plotData.empty() )
        {
            std::cout << "Nothing to plot!\n";
        }
        else
        {
            for ( auto pd : plotData)
            {
                file << " \n";
                if ( pd.pointsFill )
                {
                    for ( size_t i=0; i<pd.rs.size(); ++i )
                    {
                        file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(pd.rs[i],pd.thetas[i]))
                            << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(pd.rs[i],pd.thetas[i]))
                            << "\" r=\"" << pd.pointsRadius
                            << "\" stroke=\"" << pd.pointsBorderColor
                            << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                            << "\" stroke-width=\"" << pd.pointsBorderWidth
                            << "\" fill=\"" << pd.pointsFillColor
                            << "\" fill-opacity=\"" << pd.pointsFillOpacity
                            << "\" />\n";
                    }
                }
                else
                {
                    for ( size_t i=0; i<rVals.size(); ++i )
                    {
                        file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(pd.rs[i],pd.thetas[i]))
                            << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(pd.rs[i],pd.thetas[i]))
                            << "\" r=\"" << pd.pointsRadius
                            << "\" stroke=\"" << pd.pointsBorderColor
                            << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                            << "\" fill=\"none\" stroke-width=\"" << pd.pointsBorderWidth << "\" />\n";
                    }
                }
            }
            std::cout << "Plot generated and compiled to: " << filename << '\n';
        }
    }
    else
    {
        if ( plotData.empty() )
        {
            file << " \n";
            if ( PlotPolar::plotPointsFill )
            {
                for ( size_t i=0; i<rVals.size(); ++i )
                {
                    file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(rVals[i],thetaVals[i]))
                        << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(rVals[i],thetaVals[i]))
                        << "\" r=\"" << PlotPolar::plotPointsRadius
                        << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                        << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                        << "\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth
                        << "\" fill=\"" << PlotPolar::plotPointsFillColor
                        << "\" fill-opacity=\"" << PlotPolar::plotPointsFillOpacity
                        << "\" />\n";
                }
            }
            else
            {
                for ( size_t i=0; i<rVals.size(); ++i )
                {
                    file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(rVals[i],thetaVals[i]))
                        << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(rVals[i],thetaVals[i]))
                        << "\" r=\"" << PlotPolar::plotPointsRadius
                        << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                        << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                        << "\" fill=\"none\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth << "\" />\n";
                }
            }
            std::cout << "Plot generated and compiled to: " << filename << '\n';
        }
        else
        {
            file << " \n";
            if ( PlotPolar::plotPointsFill )
            {
                for ( size_t i=0; i<rVals.size(); ++i )
                {
                    file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(rVals[i],thetaVals[i]))
                        << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(rVals[i],thetaVals[i]))
                        << "\" r=\"" << PlotPolar::plotPointsRadius
                        << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                        << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                        << "\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth
                        << "\" fill=\"" << PlotPolar::plotPointsFillColor
                        << "\" fill-opacity=\"" << PlotPolar::plotPointsFillOpacity
                        << "\" />\n";
                }
            }
            else
            {
                for ( size_t i=0; i<rVals.size(); ++i )
                {
                    file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(rVals[i],thetaVals[i]))
                        << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(rVals[i],thetaVals[i]))
                        << "\" r=\"" << PlotPolar::plotPointsRadius
                        << "\" stroke=\"" << PlotPolar::plotPointsBorderColor
                        << "\" stroke-opacity=\"" << PlotPolar::plotPointsBorderOpacity
                        << "\" fill=\"none\" stroke-width=\"" << PlotPolar::plotPointsBorderWidth << "\" />\n";
                }
            }
            for ( auto pd : plotData)
            {
                file << " \n";
                if ( pd.pointsFill )
                {
                    for ( size_t i=0; i<pd.rs.size(); ++i )
                    {
                        file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(pd.rs[i],pd.thetas[i]))
                            << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(pd.rs[i],pd.thetas[i]))
                            << "\" r=\"" << pd.pointsRadius
                            << "\" stroke=\"" << pd.pointsBorderColor
                            << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                            << "\" stroke-width=\"" << pd.pointsBorderWidth
                            << "\" fill=\"" << pd.pointsFillColor
                            << "\" fill-opacity=\"" << pd.pointsFillOpacity
                            << "\" />\n";
                    }
                }
                else
                {
                    for ( size_t i=0; i<rVals.size(); ++i )
                    {
                        file << " <circle cx=\"" << PlotPolar::toPixelX(PlotPolar::rttox(pd.rs[i],pd.thetas[i]))
                            << "\" cy=\"" << PlotPolar::toPixelY(PlotPolar::rttoy(pd.rs[i],pd.thetas[i]))
                            << "\" r=\"" << pd.pointsRadius
                            << "\" stroke=\"" << pd.pointsBorderColor
                            << "\" stroke-opacity=\"" << pd.pointsBorderOpacity
                            << "\" fill=\"none\" stroke-width=\"" << pd.pointsBorderWidth << "\" />\n";
                    }
                }
            }
            std::cout << "Plot generated and compiled to: " << filename << '\n';
        }
    }
    // Write Title
    if ( PlotPolar::plotTitle )
    {
        writeplotTitle(file);
    }
    if ( PlotPolar::plotLegend )
    {
        std::cout << "Drawing the legend!\n";
        drawLegend(file);
    }
    // End
    file << "</svg>\n";
    file.close();
}
