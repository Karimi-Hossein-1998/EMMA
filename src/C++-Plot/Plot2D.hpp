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

inline double estimateStringWidth(const std::string& text, double fontsize)
{
    std::string text_ = text;
    size_t pos = 0;

    // --- First loop: Replace HTML entities like &#XXXX; with '^' ---
    while (true)
    {
        size_t startPos = text_.find("&#", pos);
        if (startPos == std::string::npos)
        {
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
        if (startPos == std::string::npos)
        {
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

enum class LineStyle
{
    Solid,
    Dashed,
    Dotted,
    DashDot
};

class Plot2DData
{
    public:
        dVec xs;
        dVec ys;
        LineStyle lineStyle;
        std::string lineColor;
        double lineWidth;
        double lineOpacity = 1.0;
        std::string lineLabel = "";
        bool errorPlot = false;
        dVec errors;
        dVec ysErrorsPositive;
        dVec ysErrorsNegative;

        Plot2DData(const dVec& x, const dVec& y, LineStyle ls, size_t R, size_t G, size_t B, double lw, double op = 1.0, std::string ll = "", bool ep=false, const dVec& errs={})
        : xs(x), ys(y), lineStyle(ls), lineWidth(lw), lineOpacity(op), lineLabel(ll), errorPlot(ep)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                lineColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                lineColor = "#000000";
                std::cout << "Values were out of bound! Defaulted to black for lineColor!\n";
            }
            if (errs.size() == ys.size())
            {
                errors = errs;
                for ( size_t i=0; i<ys.size(); ++i )
                {
                    ysErrorsPositive.push_back(ys[i]+errors[i]);
                    ysErrorsNegative.push_back(ys[i]-errors[i]);
                }
            }
            else if (errorPlot)
            {
				throw std::runtime_error("Errors must be the same shape and size as the y values!\n");
			}
        };
};

inline std::string dashArrayAttribute(LineStyle ls, double w)
{
    std::string dashArrayAttr;
    // writing original data
    switch (ls)
    {
        case LineStyle::Dashed:
            dashArrayAttr = " stroke-dasharray=\""+std::to_string(2.0*w)+","+std::to_string(w)+"\"";
            break;
        case LineStyle::Dotted:
            dashArrayAttr = " stroke-dasharray=\""+std::to_string(0.05*w)+","+std::to_string(w)+"\"";
            break;
        case LineStyle::DashDot:
            dashArrayAttr = " stroke-dasharray=\""+std::to_string(2.0*w)+","+std::to_string(w)
            +std::to_string(0.05*w)+","+std::to_string(w)+"\"";
            break;
        case LineStyle::Solid:
        default:
            dashArrayAttr = "";
            break;
    }
	return dashArrayAttr;
}

enum class LegendPos
{
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight
};

class Plot2D
{
    private:
        bool forcedXMinMax = false;
        bool forcedYMinMax = false;
        // Error
        dVec errors;
        bool errorPlot = false;
        dVec yValsErrorsPositive;
        dVec yValsErrorsNegative;

        double epsilon = 1e-10;
        // Set in initiatioin
        Vec<Plot2DData> plotData;
        size_t width;
        size_t height;
        dVec   xVals;
        dVec   yVals;
        size_t padding;
        double xMax;
        double xMin;
        double xRange;
        double yMax;
        double yMin;
        double yRange;
        double xScale;
        double yScale;
        // Actual Drawn width and height;
        size_t drawW;
        size_t drawH;
        bool   xAxisVisible;
        bool   yAxisVisible;
        // Set inside `drawTicks` function
        dVec xTickPixelsX;
        dVec xTickPointsX;
        dVec yTickPixelsY;
        dVec yTickPointsY;
        double gLeft;
        double gRight;
        double gBottom;
        double gTop;
        size_t numYTicks;
        double xStep;
        double yStep;
        // Set manually
        bool   border                        = true;
        bool   axis                          = true;
        bool   ticks                         = true;
        bool   grid                          = true;
        bool   axisLabel                     = true;
        bool   autoTicks                     = true;
        bool   plotTitle                     = false;
        bool   plotLegend                    = false;
        size_t numXTicks                     = 11;
        double tickMarkSize                  = 6.0;
        double plotPad                       = 4.0;
        double borderLineWidth               = 1.5;
        double axisLineWidth                 = 2.5;
        double gridLineWidth                 = 0.5;
        double plotNumeralsFontSize          = 15.0;
        double plotAxisLabelFontSize         = 16.0;
        double plotTitleFontSize             = 20.0;
        double legendFontSize                = 10;
        double plotLineWidth                 = 1.0;
        size_t plotNumeralsSignificantDigits = 3;
        LineStyle plotLineStyle = LineStyle::Solid;
        LegendPos legendPos = LegendPos::TopLeft;
        double plotLineOpacity = 1.0;
        std::string xAxisLabel = "X";
        std::string yAxisLabel = "Y";
        std::string plotTitleText = "Title!";
        // ----- Colors ----- //
        std::string backColor          = "#ffffff";
        std::string borderLineColor    = "#000000";
        std::string axisLineColor      = "#000000";
        std::string gridLineColor      = "#101010";
        std::string plotNumeralsColor  = "#000000";
        std::string plotLineColor      = "#0000ff";
        std::string plotAxisLabelColor = "#1f1f1f";
        std::string plotTitleColor     = "#000000";
        std::string plotLineLabel      = "";

        // Segmented plotting
        std::string canvasString = "";
        std::string graphString = "";
        std::string legendBoxString = "";
        // Pont -> Pixel & Pixel -> Point
        inline double toPixelX( const double x) const {return padding + (x-xMin)*xScale;};
        inline double toPixelY( const double y) const {return (height-padding) - (y-yMin)*yScale;};
        inline double toPointX( const double x) const {return (x-padding)/xScale + xMin;};
        inline double toPointY( const double y) const {return yMin - (y-height+padding)/yScale;};
        // Plot options
        inline void drawTicks(std::stringstream& file);
        inline void drawGrid(std::stringstream& file);
        inline void writeNumbers(std::stringstream& file);
        inline void writeAxisLabels(std::stringstream& file);
        inline void writeplotTitle(std::stringstream& file);
        inline void drawLegend();
        inline double calcAutoXTicksStep()
        {
            double rawStep  = xRange / (numXTicks - 1);
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
        inline double calcAutoYTicksStep()
        {
            double rawStep  = yRange / (numYTicks - 1);
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

        double yNumeralLongest = 0.0;
        inline void adjustRange()
        {
            xRange = xMax - xMin;
            yRange = yMax - yMin;
            if ( xRange<2.0*epsilon )
            {
                xMin -= epsilon;
                xMax += epsilon;
                xRange = xMax-xMin;
            }
            if ( yRange<2.0*epsilon )
            {
                yMin -= epsilon;
                yMax += epsilon;
                yRange = yMax-yMin;
            }
            // Calculating Scalings;
            xScale = (std::abs(xRange) > epsilon)? (drawW/xRange) : 1.0;
            yScale = (std::abs(yRange) > epsilon)? (drawH/yRange) : 1.0;
            xAxisVisible = ((yMax>0) && (yMin<=0)) || ((yMax>=0) && (yMin<0));
            yAxisVisible = ((xMax>0) && (xMin<=0)) || ((xMax>=0) && (xMin<0));
        };
        inline void adjustMinMax(const dVec& x, const dVec& y)
        {
            if (!forcedXMinMax)
            {
                xMin = std::min(*std::min_element(x.begin(),x.end()), xMin);
                xMax = std::max(*std::max_element(x.begin(),x.end()), xMax);
            }
            if (!forcedYMinMax)
            {
                yMin = std::min(*std::min_element(y.begin(),y.end()), yMin);
                yMax = std::max(*std::max_element(y.begin(),y.end()), yMax);
            }
        };
        inline void calculateTickPlacement();
    public:
        Plot2D(const Vec<Plot2DData>& data, size_t w = 1000, size_t h = 1000, size_t pad = 100)
        : width(w), height(h), padding(pad)
        {
            if ( data.empty() )
                throw std::runtime_error("The data is empty! You must provide at least one set of data.");
            if ( (2*padding>=width) || (2*padding>=height) )
                throw std::runtime_error("Padding is two large! It exceeds the height or width parameters.");
            for ( const auto& pd : data )
            {
                if ( (pd.xs.size() != pd.ys.size()) || (pd.xs.empty()) || (pd.ys.empty()) )
                    continue;
                if ( pd.errorPlot && (pd.errors.size()!=pd.ys.size()) )
                    continue;
				plotData.push_back(pd);
            }
        };
        Plot2D(const dVec& xs, const dVec& ys, size_t w = 1000, size_t h = 1000, size_t pad = 100)
        : width(w), height(h), xVals(xs), yVals(ys), padding(pad)
        {
            // Simple Safety Checks!
            if ( xVals.empty() )
                throw std::runtime_error("Empty data! Cannot plot empty dataset!");
            if ( xVals.size() != yVals.size() )
                throw std::runtime_error("Dimensional mismatch! 'x' and 'y' must have the same dimensions!");
            if ( (2*padding>=width) || (2*padding>=height) )
                throw std::runtime_error("Padding is two large! It exceeds the height or width parameters.");
        };
        inline void Setup();
        ~Plot2D(){};
        inline void addData(const dVec& x, const dVec& y, LineStyle ls, size_t R, size_t G, size_t B, double lw, double op = 1.0, std::string ll = "", bool ep=false, const dVec& errs={})
        {
            // Simple Safety Checks!
            if ( x.size() != y.size() || x.empty() )
                throw std::runtime_error("Dimensional mismatch! 'x' and 'y' must have the same dimensions!");
            if ( ep && (errs.size()!=y.size()) )
                throw std::runtime_error("Errors must be the same shape and size as y values!\n");
            Plot2DData pd(x,y,ls,R,G,B,lw,op,ll,ep,errs);
            Plot2D::plotData.push_back(pd);
        };
        inline void addMultipleData(const Vec<Plot2DData>& plotdata)
        {
            if ( plotdata.empty() )
                throw std::runtime_error("The data is empty! You must provide at least one set of data.");
            for ( const auto& pd : plotdata )
            {
                if ( (pd.xs.size() != pd.ys.size()) || (pd.xs.empty()) ) continue;
                if ( pd.errorPlot && (pd.errors.size()!=pd.ys.size()) ) continue;
                Plot2D::plotData.push_back(pd);
            }
        };
        inline void forceXMinMax(double xMin_, double xMax_)
        {
            xMin = xMin_;
            xMax = xMax_;
            adjustRange();
            forcedXMinMax = true;
        }
        inline void forceYMinMax(double yMin_, double yMax_)
        {
            yMin = yMin_;
            yMax = yMax_;
            adjustRange();
            forcedYMinMax = true;
        }
        // Plot
        inline void plotSVG(const std::string& filename);

        inline void Canvas();
        inline void Plot();
        inline void Plot(const Plot2DData& pd);
        inline void Legend();
        inline void Save(std::string filename);
        // Setting values
        inline void setborder(bool border_) {border=border_;};
        inline void setaxis(bool axis_) {axis=axis_;};
        inline void setticks(bool ticks_) {ticks=ticks_;};
        inline void setgrid(bool grid_) {grid=grid_;};
        inline void setplotTitle(bool title_) {plotTitle=title_;};
        inline void setnumXTicks(size_t nt) {numXTicks=nt;};
        inline void settickMarkSize(double ts) {tickMarkSize=ts;};
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
        inline void setaxisLineWidth(double alw) {axisLineWidth=alw;};
        inline void setaxisLineColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                axisLineColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to black axes lines...";
            }
        };
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
        inline void setplotNumeralsFontSize(double pnfs) {plotNumeralsFontSize=pnfs;};
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
        inline void setplotLineWidth(double plw) {plotLineWidth=plw;};
        inline void setplotLineColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                plotLineColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to blue plot lines...";
            }
        };
        inline void setplotNumeralsSignificantDigits(size_t pnsd) {plotNumeralsSignificantDigits=pnsd;};
        inline void setaxisLabel(bool al) {axisLabel=al;};
        inline void setxAxisLabel(std::string xal) {xAxisLabel=xal;};
        inline void setyAxisLabel(std::string yal) {yAxisLabel=yal;};
        inline void setplotAxisLabelFontSize(double palfs) {plotAxisLabelFontSize=palfs;};
        inline void setplotAxisLabelColor(size_t R, size_t G, size_t B)
        {
            if ( (R<256) && (G<256) && (B<256) )
            {
                plotAxisLabelColor = "#"+((R<16) ? "0"+std::format("{:x}",R) : std::format("{:x}",R))
                    +((G<16) ? "0"+std::format("{:x}",G) : std::format("{:x}",G))
                    +((B<16) ? "0"+std::format("{:x}",B) : std::format("{:x}",B));
            }
            else
            {
                std::cout << "`R`, `G`, and `B` must be integers from 0 upto 255!\nDefaulting to dark gray plot axis labels...";
            }
        };
        inline void setplotTitleFontSize(double ptfs) {plotTitleFontSize=ptfs;};
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
        inline void setautoTicks(bool at) {autoTicks=at;};
        inline void setplotLineStyle(LineStyle pls) {plotLineStyle=pls;};
        inline void setplotTitleText(std::string titletext) {plotTitleText=titletext;};
        inline void setplotLineOpacity(double op) {plotLineOpacity=op;};
        inline void setplotLegend(bool pl) {plotLegend=pl;};
        inline void setplotLineLabel(std::string pll) {plotLineLabel=pll;};
        inline void setlegendPos(LegendPos lp) {legendPos=lp;};
        inline void setlegendFontSize(double lfs) {legendFontSize=lfs;};
        inline void setplotErrors(const dVec& errs)
        {
            errors=errs;
            if ( errors.empty() || errors.size()!=yVals.size() )
            {
                throw std::runtime_error("Errors must have the same size as the data!\n");
            }
            for ( size_t i=0; i<yVals.size(); ++i )
            {
                yValsErrorsNegative.push_back(yVals[i]-errors[i]);
                yValsErrorsPositive.push_back(yVals[i]+errors[i]);
            }
        };
        inline void seterrorPlot(bool ep) {errorPlot=ep;};
};

inline void Plot2D::calculateTickPlacement()
{
    if ( Plot2D::autoTicks ) // Clean Tick Placing
    {
        Plot2D::xStep     = Plot2D::calcAutoXTicksStep();
        Plot2D::numYTicks = (Plot2D::drawW==Plot2D::drawH) ? Plot2D::numXTicks : static_cast<size_t>(Plot2D::numXTicks*Plot2D::drawH/Plot2D::drawW);
        Plot2D::yStep     = Plot2D::calcAutoYTicksStep();
        // finding first clean point of y tick.
        double firstY = std::ceil(Plot2D::yMin/Plot2D::yStep)*Plot2D::yStep;
        for ( double val=firstY; val<=Plot2D::yMax+(Plot2D::yStep*1.0e-5); val+=yStep )
        {
            Plot2D::yTickPixelsY.push_back(Plot2D::toPixelY(val));
            Plot2D::yTickPointsY.push_back(val);
        }
        // finding first clean point of x tick.
        double firstX = std::ceil(Plot2D::xMin/Plot2D::xStep)*Plot2D::xStep;
        for ( double val=firstX; val<=Plot2D::xMax+(Plot2D::xStep*1.0e-5); val+=xStep )
        {
            Plot2D::xTickPixelsX.push_back(Plot2D::toPixelX(val));
            Plot2D::xTickPointsX.push_back(val);
        }
    }
    else // Rough Tick placing
    {
        Plot2D::xStep     = (Plot2D::xRange) / (Plot2D::numXTicks - 1);
        Plot2D::numYTicks = (Plot2D::drawW==Plot2D::drawH) ? Plot2D::numXTicks : static_cast<size_t>(Plot2D::numXTicks*Plot2D::drawH/Plot2D::drawW);
        Plot2D::yStep     = (Plot2D::yRange) / (Plot2D::numYTicks-1);
        // yTicks
        if ( Plot2D::xAxisVisible )
        {
            for ( double tickyVal=0.0; tickyVal<=Plot2D::yMax; tickyVal+=Plot2D::yStep)
            {
                Plot2D::yTickPixelsY.push_back(Plot2D::toPixelY(tickyVal));
                Plot2D::yTickPointsY.push_back(tickyVal);
            }
            for ( double tickyVal=0.0; tickyVal>=Plot2D::yMin; tickyVal-=Plot2D::yStep)
            {
                Plot2D::yTickPixelsY.push_back(Plot2D::toPixelY(tickyVal));
                Plot2D::yTickPointsY.push_back(tickyVal);
            }
        }
        else
        {
            for ( size_t i=0; i<Plot2D::numYTicks; ++i)
            {
                Plot2D::yTickPixelsY.push_back(Plot2D::toPixelY(Plot2D::yMin+i*Plot2D::yStep));
                Plot2D::yTickPointsY.push_back(Plot2D::yMin+i*Plot2D::yStep);
            }
            // xTicks
            if ( Plot2D::yAxisVisible )
            {
                for ( double tickxVal=0.0; tickxVal<=Plot2D::xMax; tickxVal+=Plot2D::xStep)
                {
                    Plot2D::xTickPixelsX.push_back(Plot2D::toPixelX(tickxVal));
                    Plot2D::xTickPointsX.push_back(tickxVal);
                }
                for ( double tickxVal=0.0; tickxVal>=Plot2D::xMin; tickxVal-=Plot2D::xStep)
                {
                    Plot2D::xTickPixelsX.push_back(Plot2D::toPixelX(tickxVal));
                    Plot2D::xTickPointsX.push_back(tickxVal);
                }
            }
            else
            {
                for ( size_t i=0; i<Plot2D::numXTicks; ++i)
                {
                    Plot2D::xTickPixelsX.push_back(Plot2D::toPixelX(Plot2D::xMin+i*Plot2D::xStep));
                    Plot2D::xTickPointsX.push_back(Plot2D::xMin+i*Plot2D::xStep);
                }
            }
        }
    }
}

inline void Plot2D::Setup()
{
    // Actual Drawn width and height;
    drawW = width - (padding*2);
    drawH = height - (padding*2);
    // Setting up the canvas details
    if (!xVals.empty())
    {
        if (!forcedXMinMax)
        {
            xMax = *std::max_element(xVals.begin(), xVals.end());
            xMin = *std::min_element(xVals.begin(), xVals.end());
        }
        if (!forcedYMinMax)
        {
            yMax = *std::max_element(yVals.begin(), yVals.end());
            yMin = *std::min_element(yVals.begin(), yVals.end());
        }
    }
    if ( !plotData.empty())
    {
        for ( const auto& pd : plotData )
        {
            adjustMinMax(pd.xs,pd.ys);
        }
    }
    adjustRange();
    // Marking corners
    Plot2D::gLeft   = Plot2D::toPixelX(Plot2D::xMin);
    Plot2D::gRight  = Plot2D::toPixelX(Plot2D::xMax);
    Plot2D::gTop    = Plot2D::toPixelY(Plot2D::yMax);
    Plot2D::gBottom = Plot2D::toPixelY(Plot2D::yMin);
    calculateTickPlacement();
}

inline void Plot2D::drawTicks(std::stringstream& file)
{
    double xticksY = Plot2D::gBottom+Plot2D::plotPad;
    double yticksX = Plot2D::gLeft-Plot2D::plotPad;
    double xOrigin = Plot2D::toPixelX(0.0);
    double yOrigin = Plot2D::toPixelY(0.0);
    file << " \n";
    // yTicks
    for ( size_t i=0; i<Plot2D::numYTicks; ++i)
    {
        if (Plot2D::yAxisVisible)
        {
            file << " <line x1=\"" << xOrigin-Plot2D::tickMarkSize/2.0 << "\" y1=\"" << yTickPixelsY[i]
                << "\" x2=\"" << xOrigin+Plot2D::tickMarkSize/2.0 << "\" y2=\"" << yTickPixelsY[i]
                << "\" stroke=\"" << Plot2D::axisLineColor
                << "\" stroke-width=\"" << Plot2D::axisLineWidth << "\" stroke-linecap=\"round\" />\n";
        }
        file << " <line x1=\"" << yticksX-Plot2D::tickMarkSize << "\" y1=\"" << yTickPixelsY[i]
            << "\" x2=\"" << yticksX << "\" y2=\"" << yTickPixelsY[i]
            << "\" stroke=\"" << Plot2D::borderLineColor << "\" stroke-width=\"" << Plot2D::borderLineWidth << "\" stroke-linecap=\"round\" />\n";
    }
    // xTicks
    file << " \n";
    for ( size_t i=0; i<Plot2D::numXTicks; ++i)
    {
        if (Plot2D::xAxisVisible)
        {
            file << " <line x1=\"" << xTickPixelsX[i] << "\" y1=\"" << yOrigin-Plot2D::tickMarkSize/2.0
                << "\" x2=\"" << xTickPixelsX[i] << "\" y2=\"" << yOrigin+Plot2D::tickMarkSize/2.0
                << "\" stroke=\"" << Plot2D::axisLineColor << "\" stroke-width=\"" << Plot2D::axisLineWidth << "\" stroke-linecap=\"round\" />\n";
        }
        file << " <line x1=\"" << xTickPixelsX[i] << "\" y1=\"" << xticksY
            << "\" x2=\"" << xTickPixelsX[i] << "\" y2=\"" << xticksY+Plot2D::tickMarkSize
            << "\" stroke=\"" << Plot2D::borderLineColor << "\" stroke-width=\"" << Plot2D::borderLineWidth << "\" stroke-linecap=\"round\" />\n";
    }
}

inline void Plot2D::drawGrid(std::stringstream& file)
{
    if (!Plot2D::ticks)
    {

    }
    file << " \n";
    for ( size_t i=0; i<Plot2D::xTickPointsX.size(); ++i)
    {
        file << " <line x1=\"" << Plot2D::xTickPixelsX[i] << "\" y1=\"" << Plot2D::gTop-Plot2D::plotPad
            << "\" x2=\"" << Plot2D::xTickPixelsX[i] << "\" y2=\"" << Plot2D::gBottom+Plot2D::plotPad
            << "\" stroke=\"" << Plot2D::gridLineColor << "\" stroke-width=\"" << Plot2D::gridLineWidth << "\" stroke-linecap=\"round\" />\n";
    }
    for ( size_t i=0; i<Plot2D::yTickPointsY.size(); ++i )
    {
        file << " <line x1=\"" << Plot2D::gLeft-Plot2D::plotPad << "\" y1=\"" << Plot2D::yTickPixelsY[i]
            << "\" x2=\"" << Plot2D::gRight+Plot2D::plotPad << "\" y2=\"" << Plot2D::yTickPixelsY[i]
            << "\" stroke=\"" << Plot2D::gridLineColor << "\" stroke-width=\"" << Plot2D::gridLineWidth << "\" stroke-linecap=\"round\" />\n";
    }
}

inline void Plot2D::writeNumbers(std::stringstream& file)
{
    file << " \n";
    file << " <g font-family=\"sans-serif\" font-size=\"" << Plot2D::plotNumeralsFontSize
         << "\" fill=\"" << Plot2D::plotNumeralsColor << "\">\n";

    // --- X-AXIS NUMERALS ---
    for ( size_t i=0; i<Plot2D::xTickPointsX.size(); ++i)
    {
        std::stringstream ss;
        // defaultfloat strips trailing zeros and manages scientific boundaries dynamically
        // 3 significant digits means: 1234, 123.4, 12.34, 1.234, or 1.234e+05
        ss << std::defaultfloat << std::setprecision(Plot2D::plotNumeralsSignificantDigits) << Plot2D::xTickPointsX[i];
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

        file << "   <text x=\"" << Plot2D::xTickPixelsX[i]
             << "\" y=\"" << Plot2D::gBottom+Plot2D::tickMarkSize+Plot2D::plotPad+Plot2D::plotNumeralsFontSize*0.5
             << "\" text-anchor=\"middle\" dominant-baseline=\"central\">"
             << labelStr << "</text>\n";
    }

    // --- Y-AXIS NUMERALS ---
    for ( size_t i=0; i<Plot2D::yTickPointsY.size(); ++i )
    {
        std::stringstream ss;
        ss << std::defaultfloat << std::setprecision(Plot2D::plotNumeralsSignificantDigits) << Plot2D::yTickPointsY[i];
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

        Plot2D::yNumeralLongest = std::max(estimateStringWidth(labelStr,Plot2D::plotNumeralsFontSize),yNumeralLongest);
        // Added dominant-baseline="central" to keep numbers aligned with the tick line rows
        file << "   <text x=\"" << Plot2D::gLeft-Plot2D::tickMarkSize-Plot2D::plotPad-Plot2D::plotNumeralsFontSize*0.25
             << "\" y=\"" << Plot2D::yTickPixelsY[i]
             << "\" text-anchor=\"end\" dominant-baseline=\"central\">"
             << labelStr << "</text>\n";
    }
    file << " </g>\n";
}

inline void Plot2D::writeAxisLabels(std::stringstream& file)
{
    // --- AXIS TITLES / LABELS ---
    // file << " \n";
    file << " <g font-family=\"sans-serif\" font-size=\"" << Plot2D::plotAxisLabelFontSize
         << "\" fill=\"" << Plot2D::plotAxisLabelColor << "\">\n";

    // 1. X-Axis Label: Centered horizontally beneath the numbers row
    double xAxisLabelX = Plot2D::width / 2.0;
    // Push it down further past the numerals padding
    double xAxisLabelY = Plot2D::gBottom+Plot2D::tickMarkSize+Plot2D::plotPad
    					+Plot2D::plotNumeralsFontSize+Plot2D::plotAxisLabelFontSize*0.5;
    file << "   \n";
    file << "   <text x=\"" << xAxisLabelX << "\" y=\"" << xAxisLabelY
        << "\" text-anchor=\"middle\" dominant-baseline=\"central\" font-weight=\"bold\">"
        << Plot2D::xAxisLabel << "</text>\n";

    // 2. Y-Axis Label: Placed to the left of the Y-axis numbers and rotated -90 degrees
    // Calculate a safe padding slot out to the left of the border box
    double yAxisLabelX = 0.0;
    if (yNumeralLongest==0)
    {
        yAxisLabelX = Plot2D::gLeft - Plot2D::tickMarkSize
                    - (1.1+Plot2D::plotNumeralsSignificantDigits)*Plot2D::plotNumeralsFontSize
                    - Plot2D::plotAxisLabelFontSize*0.5;
    }
    else
    {
        yAxisLabelX = Plot2D::gLeft - Plot2D::tickMarkSize - Plot2D::plotPad - Plot2D::plotNumeralsFontSize*0.3
                    - Plot2D::yNumeralLongest
                    - Plot2D::plotAxisLabelFontSize*0.5;
    }
    // double yAxisLabelX = Plot2D::gLeft - 4.0*(Plot2D::tickMarkSize + Plot2D::plotPad) + 2.0*Plot2D::plotNumeralsFontSize;
    double yAxisLabelY = Plot2D::height / 2.0; // Vertically centered on the plot wall

    file << "   \n";
    // We pivot exactly around (yLabelX, yLabelY) to keep the text perfectly positioned
    file << "   <text x=\"" << yAxisLabelX << "\" y=\"" << yAxisLabelY
        << "\" transform=\"rotate(-90, " << yAxisLabelX << ", " << yAxisLabelY << ")"
        << "\" text-anchor=\"middle\" dominant-baseline=\"central\" font-weight=\"bold\">"
        << Plot2D::yAxisLabel << "</text>\n";
    file << " </g>\n";
}

inline void Plot2D::writeplotTitle(std::stringstream& file)
{
    file << " <g font-family=\"sans-serif\" font-size=\"" << Plot2D::plotTitleFontSize
         << "\" fill=\"" << Plot2D::plotTitleColor << "\">\n";
    file << "   <text x=\"" << Plot2D::width / 2.0 << "\" y=\"" << Plot2D::plotTitleFontSize+10.0
        << "\" text-anchor=\"middle\" dominant-baseline=\"central\" font-weight=\"bold\">"
        << Plot2D::plotTitleText << "</text>\n";
    file << " </g>\n";
}

inline void Plot2D::Canvas()
{
    std::stringstream file;

    // Background
    file << "<svg width=\"" << Plot2D::width << "\" height=\"" << Plot2D::height
        << "\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:" << Plot2D::backColor <<";\">\n";

    file << " \n";
    file << " \n";
    // Painting background
    file << " <rect width=\"" << Plot2D::width << "\" height=\"" << Plot2D::height <<
    	"\" fill=\"" << Plot2D::backColor << "\" />\n";
    // Drawing Borders
    if ( border )
    {
        file << " \n";
        file << " <rect x=\"" << Plot2D::gLeft-Plot2D::plotPad << "\" y=\"" << Plot2D::gTop-Plot2D::plotPad
            << "\" width=\"" << (Plot2D::gRight-Plot2D::gLeft+2*Plot2D::plotPad)
            << "\" height=\"" << (Plot2D::gBottom-Plot2D::gTop+2*Plot2D::plotPad)
            << "\" rx=\"" << Plot2D::plotPad << "\" ry=\"" << Plot2D::plotPad
            << "\" fill=\"none\" stroke=\"" << Plot2D::borderLineColor << "\" stroke-width=\"" << Plot2D::borderLineWidth << "\" />\n";
    }
    // Drawing axes lines;
    if ( axis )
    {
        // Drawing xAxis
        // Checking origin
        file << " \n";
        if ( Plot2D::xAxisVisible )
        {
            double xAxisX1 = Plot2D::gLeft-Plot2D::plotPad;
            double xAxisX2 = Plot2D::gRight+Plot2D::plotPad;
            double xAxisY1 = Plot2D::toPixelY(0.0);
            double xAxisY2 = Plot2D::toPixelY(0.0);
            file << " <line x1=\"" << xAxisX1 << "\" y1=\"" << xAxisY1
                << "\" x2=\"" << xAxisX2 << "\" y2=\"" << xAxisY2
                << "\" stroke=\"" << Plot2D::axisLineColor << "\" stroke-width=\"" << Plot2D::axisLineWidth << "\" stroke-linecap=\"round\" />\n";
        }
        // Drawing yAxis
        // Checking origin
        if ( Plot2D::yAxisVisible )
        {
            double yAxisX1 = Plot2D::toPixelX(0.0);
            double yAxisX2 = Plot2D::toPixelX(0.0);
            double yAxisY1 = gBottom+Plot2D::plotPad;
            double yAxisY2 = gTop-Plot2D::plotPad;
            file << " <line x1=\"" << yAxisX1 << "\" y1=\"" << yAxisY1
                << "\" x2=\"" << yAxisX2 << "\" y2=\"" << yAxisY2
                << "\" stroke=\"" << Plot2D::axisLineColor << "\" stroke-width=\"" << Plot2D::axisLineWidth << "\" stroke-linecap=\"round\" />\n";
        }
    }

    // Drawing TickMarks
    if ( Plot2D::ticks )
    {
        Plot2D::drawTicks(file);
    }
    // Writing tickNums
    // Plot2D::writeNumbers(file);
    if ( Plot2D::grid )
    {
        Plot2D::drawGrid(file);
    }
    Plot2D::writeNumbers(file);
    // Adding a clipping mechanism
    file << " \n";
    file << " <defs>\n  <clipPath id=\"rect-plot-clip\">\n   <rect x=\"" << gLeft
    	<< "\" y=\"" << gTop << "\" width=\"" << drawW << "\" height=\"" << drawH << "\" />\n  </clipPath>\n </defs>\n";
    // Dummy rectangle to force bounding box
    file << " <rect x=\"" << gLeft-1 << "\" y=\"" << gTop-1 << "\" width=\"" << drawW+2 << "\" height=\"" << drawH+2
	    << "\" fill=\"none\" pointer-events=\"none\" clip-path=\"url(rect-plot-clip)\" />";
    // Writing string to the class.
    Plot2D::canvasString = file.str();
}

inline void Plot2D::Plot()
{
    std::stringstream file;
    double xp=0.0,yp=0.0,xpPre=0.0,ypPre=0.0;
    double yContinuityThreshold = drawH*0.5;
    // Plotting the xVals, yVals
    if ( !xVals.empty() )
    {
        file << "\n";
        // The data points.
        std::string dashArrayAttr = dashArrayAttribute(Plot::plotLineStyle, Plot2D::plotLineWidth);
        std::stringstream errorData;
        std::stringstream errorLinePositive;
        std::stringstream errorLineNegative;
        if (errorPlot)
        {
            xp = std::clamp(Plot2D::toPixelX(xVals[0]),-1.0,width+1.0);
            yp = std::clamp(Plot2D::toPixelY(yValsErrorsPositive[0]),-1.0,height+1.0);
            errorData << "M " << xp << " " << yp;
            errorLinePositive << "M " << xp << " " << yp;
            yp = std::clamp(Plot2D::toPixelY(yValsErrorsNegative[0]),-1.0,height+1.0);
            errorLineNegative << "M " << xp << " " << yp;
            for ( size_t i=1; i<Plot2D::xVals.size(); ++i )
            {
                xp = std::clamp(Plot2D::toPixelX(xVals[i]),-1.0,width+1.0);
                yp = std::clamp(Plot2D::toPixelY(yValsErrorsPositive[i]),-1.0,height+1.0);
                errorData << " L " << xp << " " << yp;
                errorLinePositive << " L " << xp << " " << yp;
                yp = std::clamp(Plot2D::toPixelY(yValsErrorsNegative[i]),-1.0,height+1.0);
                errorLineNegative << " L " << xp << " " << yp;
            }
            for ( size_t i=Plot2D::xVals.size()-1; i>0; --i )
            {
                xp = std::clamp(Plot2D::toPixelX(xVals[i]),-1.0,width+1.0);
                yp = std::clamp(Plot2D::toPixelY(yValsErrorsNegative[i]),-1.0,height+1.0);
                errorData << " L " << xp << " " << yp;
            }
            errorData << " Z";
            file << " <path d=\"" << errorData.str()
                << "\" fill=\"" << Plot2D::plotLineColor
                << "\" fill-opacity=\"" << Plot2D::plotLineOpacity*0.4
                << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
            file << " <path d=\"" << errorLinePositive.str() << "\" fill=\"none"
                << "\" stroke=\"" << Plot2D::plotLineColor
                << "\" stroke-opacity=\"" << Plot2D::plotLineOpacity*0.5
                << "\" stroke-width=\"" << Plot2D::plotLineWidth*0.5
                << "\" stroke-linecap=\"round"
                << "\" stroke-linejoin=\"round\""
                << dashArrayAttr
                << " clip-path=\"url(#rect-plot-clip)\" />\n";
            file << " <path d=\"" << errorLineNegative.str() << "\" fill=\"none"
                << "\" stroke=\"" << Plot2D::plotLineColor
                << "\" stroke-opacity=\"" << Plot2D::plotLineOpacity*0.5
                << "\" stroke-width=\"" << Plot2D::plotLineWidth*0.5
                << "\" stroke-linecap=\"round"
                << "\" stroke-linejoin=\"round\""
                << dashArrayAttr
                << " clip-path=\"url(#rect-plot-clip)\" />\n";
        }
        if (xVals.size()>=std::floor(Plot2D::drawW/(13*Plot2D::plotLineWidth)))
        {
            std::stringstream pathData;
            xp = std::clamp(Plot2D::toPixelX(xVals[0]),-1.0,width+1.0);
            yp = std::clamp(Plot2D::toPixelY(yVals[0]),-1.0,height+1.0);
            pathData << "M " << xp << " " << yp;
            for ( size_t i=1; i<Plot2D::xVals.size(); ++i )
            {
                xpPre = xp;
                ypPre = yp;
                xp = std::clamp(Plot2D::toPixelX(xVals[i]),-1.0,width+1.0);
                yp = std::clamp(Plot2D::toPixelY(yVals[i]),-1.0,height+1.0);
                if (std::abs(ypPre-yp)>=yContinuityThreshold)
                    pathData << " M " << xp << " " << yp;
                else
                    pathData << " L " << xp << " " << yp;
            }
            file << " <path d=\"" << pathData.str() << "\" fill=\"none"
                << "\" stroke=\"" << Plot2D::plotLineColor
                << "\" stroke-opacity=\"" << Plot2D::plotLineOpacity
                << "\" stroke-width=\"" << Plot2D::plotLineWidth
                << "\" stroke-linecap=\"round"
                << "\" stroke-linejoin=\"round\""
                << dashArrayAttr << " clip-path=\"url(#rect-plot-clip)\" />\n";
        }
        else
        {
            std::stringstream pathData;
            xp = std::clamp(Plot2D::toPixelX(xVals[0]),-1.0,width+1.0);
            yp = std::clamp(Plot2D::toPixelY(yVals[0]),-1.0,height+1.0);
            pathData << "M " << xp << " " << yp;
            file << " <circle cx=\"" << xp << "\" cy=\"" << yp
                << "\" r=\"" << 1.25*Plot2D::plotLineWidth
                << "\" fill=\"" << Plot2D::plotLineColor
                << "\" fill-opacity=\"" << Plot2D::plotLineOpacity << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
            for ( size_t i=1; i<Plot2D::xVals.size(); ++i )
            {
                xpPre = xp;
                ypPre = yp;
                xp = std::clamp(Plot2D::toPixelX(xVals[i]),-1.0,width+1.0);
                yp = std::clamp(Plot2D::toPixelY(yVals[i]),-1.0,height+1.0);
                file << " <circle cx=\"" << xp << "\" cy=\"" << yp
                    << "\" r=\"" << 1.25*Plot2D::plotLineWidth
                    << "\" fill=\"" << Plot2D::plotLineColor
                    << "\" fill-opacity=\"" << Plot2D::plotLineOpacity << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
                if (std::abs(ypPre-yp)>=yContinuityThreshold)
                    pathData << " M " << xp << " " << yp;
                else
                    pathData << " L " << xp << " " << yp;
            }
            file << " <path d=\"" << pathData.str() << "\" fill=\"none"
                << "\" stroke=\"" << Plot2D::plotLineColor
                << "\" stroke-opacity=\"" << Plot2D::plotLineOpacity
                << "\" stroke-width=\"" << Plot2D::plotLineWidth
                << "\" stroke-linecap=\"round"
                << "\" stroke-linejoin=\"round\""
                << dashArrayAttr << " clip-path=\"url(#rect-plot-clip)\" />\n";
        }
    }
    // Plotting the plotData
    if ( !Plot2D::plotData.empty() )
    {
        for ( const auto& pd : plotData )
        {
            file << "\n";
            // The data points.
            std::string dashArrayAttr = dashArrayAttribute(pd.lineStyle, pd.lineWidth);
            std::stringstream errorData;
            std::stringstream errorLinePositive;
            std::stringstream errorLineNegative;
            if (pd.errorPlot)
            {
                xp = std::clamp(Plot2D::toPixelX(pd.xs[0]),-1.0,width+1.0);
                yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsPositive[0]),-1.0,height+1.0);
                errorData << "M " << xp << " " << yp;
                errorLinePositive << "M " << xp << " " << yp;
                yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsNegative[0]),-1.0,height+1.0);
                errorLineNegative << "M " << xp << " " << yp;
                for ( size_t i=1; i<pd.xs.size(); ++i )
                {
                    xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
                    yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsPositive[i]),-1.0,height+1.0);
                    errorData << " L " << xp << " " << yp;
                    errorLinePositive << " L " << xp << " " << yp;
                    yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsNegative[i]),-1.0,height+1.0);
                    errorLineNegative << " L " << xp << " " << yp;
                }
                for ( size_t i=Plot2D::xVals.size()-1; i>0; --i )
                {
                    xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
                    yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsNegative[i]),-1.0,height+1.0);
                    errorData << " L " << xp << " " << yp;
                }
                errorData << " Z";
                file << " <path d=\"" << errorData.str()
                    << "\" fill=\"" << pd.lineColor
                    << "\" fill-opacity=\"" << pd.lineOpacity*0.4
                    << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
                file << " <path d=\"" << errorLinePositive.str() << "\" fill=\"none"
                    << "\" stroke=\"" << pd.lineColor
                    << "\" stroke-opacity=\"" << pd.lineOpacity*0.5
                    << "\" stroke-width=\"" << pd.lineWidth*0.5
                    << "\" stroke-linecap=\"round"
                    << "\" stroke-linejoin=\"round\""
                    << dashArrayAttr
                    << " clip-path=\"url(#rect-plot-clip)\" />\n";
                file << " <path d=\"" << errorLineNegative.str() << "\" fill=\"none"
                    << "\" stroke=\"" << pd.lineColor
                    << "\" stroke-opacity=\"" << pd.lineOpacity*0.5
                    << "\" stroke-width=\"" << pd.lineWidth*0.5
                    << "\" stroke-linecap=\"round"
                    << "\" stroke-linejoin=\"round\""
                    << dashArrayAttr
                    << " clip-path=\"url(#rect-plot-clip)\" />\n";
            }
            if (pd.xs.size()>=std::floor(Plot2D::drawW/(13*pd.lineWidth)))
            {
                std::stringstream pathData;
                xp = std::clamp(Plot2D::toPixelX(pd.xs[0]),-1.0,width+1.0);
                yp = std::clamp(Plot2D::toPixelY(pd.ys[0]),-1.0,height+1.0);
                pathData << "M " << xp << " " << yp;
                for ( size_t i=1; i<pd.xs.size(); ++i )
                {
                    xpPre = xp;
                    ypPre = yp;
                    xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
                    yp = std::clamp(Plot2D::toPixelY(pd.ys[i]),-1.0,height+1.0);
                    if (std::abs(ypPre-yp)>=yContinuityThreshold)
                        pathData << " M " << xp << " " << yp;
                    else
                        pathData << " L " << xp << " " << yp;
                }
                file << " <path d=\"" << pathData.str() << "\" fill=\"none"
                    << "\" stroke=\"" << pd.lineColor
                    << "\" stroke-opacity=\"" << pd.lineOpacity
                    << "\" stroke-width=\"" << pd.lineWidth
                    << "\" stroke-linecap=\"round"
                    << "\" stroke-linejoin=\"round\""
                    << dashArrayAttr << " clip-path=\"url(#rect-plot-clip)\" />\n";
            }
            else
            {
                std::stringstream pathData;
                xp = std::clamp(Plot2D::toPixelX(pd.xs[0]),-1.0,width+1.0);
                yp = std::clamp(Plot2D::toPixelY(pd.ys[0]),-1.0,height+1.0);
                pathData << "M " << xp << " " << yp;
                file << " <circle cx=\"" << xp << "\" cy=\"" << yp
                    << "\" r=\"" << 1.25*pd.lineWidth
                    << "\" fill=\"" << pd.lineColor
                    << "\" fill-opacity=\"" << pd.lineOpacity << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
                for ( size_t i=1; i<pd.xs.size(); ++i )
                {
                    xpPre = xp;
                    ypPre = yp;
                    xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
                    yp = std::clamp(Plot2D::toPixelY(pd.ys[i]),-1.0,height+1.0);
                    file << " <circle cx=\"" << xp << "\" cy=\"" << yp
                        << "\" r=\"" << 1.25*pd.lineWidth
                        << "\" fill=\"" << pd.lineColor
                        << "\" fill-opacity=\"" << pd.lineOpacity << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
                    if (std::abs(ypPre-yp)>=yContinuityThreshold)
                        pathData << " M " << xp << " " << yp;
                    else
                        pathData << " L " << xp << " " << yp;
                }
                file << " <path d=\"" << pathData.str() << "\" fill=\"none"
                    << "\" stroke=\"" << pd.lineColor
                    << "\" stroke-opacity=\"" << pd.lineOpacity
                    << "\" stroke-width=\"" << pd.lineWidth
                    << "\" stroke-linecap=\"round"
                    << "\" stroke-linejoin=\"round\""
                    << dashArrayAttr << " clip-path=\"url(#rect-plot-clip)\" />\n";
            }
        }
    }
    graphString = file.str();
}

inline void Plot2D::Plot(const Plot2DData& pd)
{
    std::stringstream file;
    double xp=0.0,yp=0.0,xpPre=0.0,ypPre=0.0;
    double yContinuityThreshold = drawH*0.5;
    file << "\n";
    // Simple Safety Checks!
    if ( pd.xs.size() != pd.ys.size() || pd.xs.empty() )
        throw std::runtime_error("Dimensional mismatch! 'x' and 'y' must have the same dimensions!");
    if ( pd.errorPlot && (pd.errors.size()!=pd.ys.size()) )
        throw std::runtime_error("Errors must be the same shape and size as y values!\n");
	Plot2D::plotData.push_back(pd);
    // The data points.
    std::string dashArrayAttr = dashArrayAttribute(pd.lineStyle, pd.lineWidth);
    std::stringstream errorData;
    std::stringstream errorLinePositive;
    std::stringstream errorLineNegative;
    if (pd.errorPlot)
    {
        xp = std::clamp(Plot2D::toPixelX(pd.xs[0]),-1.0,width+1.0);
        yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsPositive[0]),-1.0,height+1.0);
        errorData << "M " << xp << " " << yp;
        errorLinePositive << "M " << xp << " " << yp;
        yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsNegative[0]),-1.0,height+1.0);
        errorLineNegative << "M " << xp << " " << yp;
        for ( size_t i=1; i<pd.xs.size(); ++i )
        {
            xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
            yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsPositive[i]),-1.0,height+1.0);
            errorData << " L " << xp << " " << yp;
            errorLinePositive << " L " << xp << " " << yp;
            yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsNegative[i]),-1.0,height+1.0);
            errorLineNegative << " L " << xp << " " << yp;
        }
        for ( size_t i=Plot2D::xVals.size()-1; i>0; --i )
        {
            xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
            yp = std::clamp(Plot2D::toPixelY(pd.ysErrorsNegative[i]),-1.0,height+1.0);
            errorData << " L " << xp << " " << yp;
        }
        errorData << " Z";
        file << " <path d=\"" << errorData.str()
            << "\" fill=\"" << pd.lineColor
            << "\" fill-opacity=\"" << pd.lineOpacity*0.4
            << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
        file << " <path d=\"" << errorLinePositive.str() << "\" fill=\"none"
            << "\" stroke=\"" << pd.lineColor
            << "\" stroke-opacity=\"" << pd.lineOpacity*0.5
            << "\" stroke-width=\"" << pd.lineWidth*0.5
            << "\" stroke-linecap=\"round"
            << "\" stroke-linejoin=\"round\""
            << dashArrayAttr
            << " clip-path=\"url(#rect-plot-clip)\" />\n";
        file << " <path d=\"" << errorLineNegative.str() << "\" fill=\"none"
            << "\" stroke=\"" << pd.lineColor
            << "\" stroke-opacity=\"" << pd.lineOpacity*0.5
            << "\" stroke-width=\"" << pd.lineWidth*0.5
            << "\" stroke-linecap=\"round"
            << "\" stroke-linejoin=\"round\""
            << dashArrayAttr
            << " clip-path=\"url(#rect-plot-clip)\" />\n";
    }
    if (pd.xs.size()>=std::floor(Plot2D::drawW/(13*pd.lineWidth)))
    {
        std::stringstream pathData;
        xp = std::clamp(Plot2D::toPixelX(pd.xs[0]),-1.0,width+1.0);
        yp = std::clamp(Plot2D::toPixelY(pd.ys[0]),-1.0,height+1.0);
        pathData << "M " << xp << " " << yp;
        for ( size_t i=1; i<pd.xs.size(); ++i )
        {
            xpPre = xp;
            ypPre = yp;
            xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
            yp = std::clamp(Plot2D::toPixelY(pd.ys[i]),-1.0,height+1.0);
            if (std::abs(ypPre-yp)>=yContinuityThreshold)
                pathData << " M " << xp << " " << yp;
            else
                pathData << " L " << xp << " " << yp;
        }
        file << " <path d=\"" << pathData.str() << "\" fill=\"none"
            << "\" stroke=\"" << pd.lineColor
            << "\" stroke-opacity=\"" << pd.lineOpacity
            << "\" stroke-width=\"" << pd.lineWidth
            << "\" stroke-linecap=\"round"
            << "\" stroke-linejoin=\"round\""
            << dashArrayAttr << " clip-path=\"url(#rect-plot-clip)\" />\n";
    }
    else
    {
        std::stringstream pathData;
        xp = std::clamp(Plot2D::toPixelX(pd.xs[0]),-1.0,width+1.0);
        yp = std::clamp(Plot2D::toPixelY(pd.ys[0]),-1.0,height+1.0);
        pathData << "M " << xp << " " << yp;
        file << " <circle cx=\"" << xp << "\" cy=\"" << yp
            << "\" r=\"" << 1.25*pd.lineWidth
            << "\" fill=\"" << pd.lineColor
            << "\" fill-opacity=\"" << pd.lineOpacity << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
        for ( size_t i=1; i<pd.xs.size(); ++i )
        {
            xpPre = xp;
            ypPre = yp;
            xp = std::clamp(Plot2D::toPixelX(pd.xs[i]),-1.0,width+1.0);
            yp = std::clamp(Plot2D::toPixelY(pd.ys[i]),-1.0,height+1.0);
            file << " <circle cx=\"" << xp << "\" cy=\"" << yp
                << "\" r=\"" << 1.25*pd.lineWidth
                << "\" fill=\"" << pd.lineColor
                << "\" fill-opacity=\"" << pd.lineOpacity << "\" clip-path=\"url(#rect-plot-clip)\" />\n";
            if (std::abs(ypPre-yp)>=yContinuityThreshold)
                pathData << " M " << xp << " " << yp;
            else
                pathData << " L " << xp << " " << yp;
        }
        file << " <path d=\"" << pathData.str() << "\" fill=\"none"
            << "\" stroke=\"" << pd.lineColor
            << "\" stroke-opacity=\"" << pd.lineOpacity
            << "\" stroke-width=\"" << pd.lineWidth
            << "\" stroke-linecap=\"round"
            << "\" stroke-linejoin=\"round\""
            << dashArrayAttr << " clip-path=\"url(#rect-plot-clip)\" />\n";
    }
    Plot2D::graphString += file.str();
}

inline void Plot2D::Legend()
{
    if (plotLegend)
    {
        drawLegend();
    }
    else
    {
        Plot2D::legendBoxString = "\n";
    }
    std::stringstream file;
    if ( Plot2D::axisLabel )
    {
        Plot2D::writeAxisLabels(file);
    }
    // Writing Title
    if ( Plot2D::plotTitle )
    {
        Plot2D::writeplotTitle(file);
    }
    legendBoxString += file.str();
}

inline void Plot2D::Save(std::string filename)
{
    std::ofstream file(filename, std::ios::trunc);
    file << canvasString;
    file << graphString;
    file << legendBoxString;
    file << "</svg>\n";
    file.close();
}

inline void Plot2D::drawLegend()
{
    std::stringstream file;
    double legendPad = 20.0;
    double legendTopPad = 20.0;
    double legendBottomPad = 20.0;
    double legendLeftPad = 20.0;
    double legendRightPad = 20.0;
    double legendWidth = 2.5*legendPad+50.0;
    double legendHeight = 2.0*legendPad;
    if ( !Plot2D::xVals.empty() )
    {
        legendWidth = std::max(estimateStringWidth(Plot2D::plotLineLabel, Plot2D::legendFontSize)
                                +2.5*legendPad+50.0,legendWidth);
        legendHeight += Plot2D::legendFontSize*2.0;
    }
    if ( !Plot2D::plotData.empty() )
    {
		for ( const auto& pd : Plot2D::plotData )
        {
            legendWidth = std::max(estimateStringWidth(pd.lineLabel, Plot2D::legendFontSize)
            						+2.5*legendPad+50.0,legendWidth);
        }
        legendHeight += Plot2D::plotData.size()*Plot2D::legendFontSize*2.0;
    }

    double legendX = 0.0;
    double legendY = 0.0;
    switch ( Plot2D::legendPos )
    {
        case LegendPos::TopLeft:
            legendX = Plot2D::gLeft+legendLeftPad;
            legendY = Plot2D::gTop+legendTopPad;
            break;
        case LegendPos::TopRight:
            legendX = Plot2D::gRight-legendRightPad-legendWidth;
            legendY = Plot2D::gTop+legendTopPad;
            break;
        case LegendPos::BottomLeft:
            legendX = Plot2D::gLeft+legendLeftPad;
            legendY = Plot2D::gBottom-legendBottomPad-legendHeight;
            break;
        case LegendPos::BottomRight:
            legendX = Plot2D::gRight-legendRightPad-legendWidth;
            legendY = Plot2D::gBottom-legendBottomPad-legendHeight;
            break;
        default:
            legendX = Plot2D::gLeft+legendLeftPad;
            legendY = Plot2D::gTop+legendTopPad;
            break;
    }
    // Drawing the empty legend box
    file << "   <rect x=\"" << legendX << "\" y=\"" << legendY
        << "\" width=\"" << legendWidth << "\" height=\"" << legendHeight
        << "\" rx=\"" << legendPad << "\" ry=\"" << legendPad
        << "\" fill=\"" << Plot2D::backColor << "\" fill-opacity=\"0.75\" stroke=\"" << Plot2D::borderLineColor
        << "\" stroke-width=\"1\" />\n";
    // size_t counter = 0;
    double itemY = 0.0;
    std::string dashArrayAttr = "";
    if ( !Plot2D::xVals.empty() )
    {
        // writing original data
        dashArrayAttr = dashArrayAttribute(Plot2D::plotLineStyle, Plot2D::plotLineWidth);
        itemY = legendY + legendPad;
        file << " <line x1=\"" << legendX+legendPad << "\" y1=\"" << itemY
            << "\" x2=\"" << legendX+legendPad+50.0 << "\" y2=\"" << itemY+2.0*Plot2D::legendFontSize
            << "\" stroke=\"" << Plot2D::plotLineColor
            << "\" stroke-opacity=\"" << Plot2D::plotLineOpacity
            << "\" stroke-width=\"" << Plot2D::plotLineWidth
            << "\" stroke-linecap=\"round"
            << "\" stroke-linejoin=\"round\""
            << dashArrayAttr << " />\n";
        file << "   <text x=\"" << legendX+50.0+1.5*legendPad << "\" y=\"" << itemY+Plot2D::legendFontSize
            << "\" font-family=\"sans-serif\" font-size=\"" << Plot2D::legendFontSize
            << "\" fill=\"" << Plot2D::plotLineColor
            << "\" text-anchor=\"start\" dominant-baseline=\"central\">"
            << Plot2D::plotLineLabel << "</text>\n";
        // counter += 1;
    }
    if ( !Plot2D::plotData.empty() )
    {
        for ( const auto& pd : Plot2D::plotData )
        {
            // itemY = legendY + counter*Plot2D::legendFontSize*2.0 + legendPad;
            itemY += Plot2D::legendFontSize*2.0;
            dashArrayAttr = dashArrayAttribute(pd.lineStyle, pd.lineWidth);
            file << " <line x1=\"" << legendX+legendPad << "\" y1=\"" << itemY
                << "\" x2=\"" << legendX+legendPad+50.0 << "\" y2=\"" << itemY+2.0*Plot2D::legendFontSize
                << "\" stroke=\"" << pd.lineColor
                << "\" stroke-opacity=\"" << pd.lineOpacity
                << "\" stroke-width=\"" << pd.lineWidth
                << "\" stroke-linecap=\"round"
                << "\" stroke-linejoin=\"round\""
                << dashArrayAttr << " />\n";
            file << "   <text x=\"" << legendX+50.0+1.5*legendPad << "\" y=\"" << itemY+Plot2D::legendFontSize
                << "\" font-family=\"sans-serif\" font-size=\"" << Plot2D::legendFontSize
                << "\" fill=\"" << pd.lineColor
                << "\" text-anchor=\"start\" dominant-baseline=\"central\">"
                << pd.lineLabel << "</text>\n";
        }
    }
    legendBoxString = file.str();
}

inline void Plot2D::plotSVG(const std::string& filename)
{
    Plot2D::Setup();
    std::ofstream file(filename, std::ios::trunc);
    Plot2D::Canvas();
    file << Plot2D::canvasString;

    Plot2D::Plot();
    file << Plot2D::graphString;

    Plot2D::Legend();
    file << Plot2D::legendBoxString;

    file << "</svg>\n";
    file.close();
}
