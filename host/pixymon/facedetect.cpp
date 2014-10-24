//
// begin license header
//
// This file is part of Pixy CMUcam5 or "Pixy" for short
//
// All Pixy source code is provided under the terms of the
// GNU General Public License v2 (http://www.gnu.org/licenses/gpl-2.0.html).
// Those wishing to use Pixy source code, software and/or
// technologies under different licensing terms should contact us at
// cmucam@cs.cmu.edu. Such licensing terms are available for
// all portions of the Pixy codebase presented here.
//
// end license header
//

#include "facedetect.h"
#include <math.h>
#include <iostream>


IntegralImage::IntegralImage(const QImage &image) :
    m_height(image.height()), //TODO change to correct height
    m_width(image.width()) //TODO change to correct width
{

    // TODO resize

    // create integral image
    // initialize to n+1xm+1 array of 32-bit ints
    for (unsigned int row = 1; row < image.height(); ++row)
    {
        for (unsigned int col = 1; col < image.height(); ++col)
        {
            // get grayscale value of pixel
            uint8_t intensity = qGray(image.pixel(col, row));

            /*
            uint32_t rowSum = m_data[(row-1)*m_width + col];
            uint32_t colSum = m_data[row*m_width + (col-1)];
            uint32_t diagSum = m_data[(row-1)*m_width + (col-1)];
            */
            uint32_t rowSum = m_data[row-1][col];
            uint32_t colSum = m_data[row][col-1];
            uint32_t diagSum = m_data[row-1][col-1];
            m_data[row][col] = intensity + rowSum + colSum - diagSum;

        }
    }
}


IntegralImage::~IntegralImage()
{
}


const uint32_t IntegralImage::operator()(const uint16_t &row, const uint16_t &col) const
{
    return m_data[row][col];
}


uint32_t IntegralImage::operator()(const uint16_t &row, const uint16_t &col)
{
    return m_data[row][col];
}


uint16_t IntegralImage::width()
{
    return m_width;
}


uint16_t IntegralImage::height()
{
    return m_height;
}


detectionLocation::detectionLocation(const uint16_t &x, const uint16_t &y, const uint16_t &w, const uint16_t &h) :
    locationX(x),
    locationY(y),
    width(w),
    height(h)
{
}


LBPFeature::LBPFeature(tinyxml2::XMLElement *featureElement, tinyxml2::XMLElement *rectanglesElement)
{
    const char * internalNodeStr = featureElement->FirstChildElement("internalNodes")->GetText();
    const char * leafValStr = featureElement->FirstChildElement("leafValues")->GetText();   

    std::vector<uint32_t> internalNodes = textToUnsignedList(internalNodeStr);

    //XXX (debug)
    std::cerr << "\t\tInternal Nodes: " << std::endl;
    unsigned int nodeIdx = 0;

    for (std::vector<uint32_t>::iterator iter = internalNodes.begin(); iter != internalNodes.end(); ++iter)
    {
        std::cerr << "\t\tNode " << nodeIdx++ << ": " << *iter << std::endl;
    }

    std::vector<double> leafValues = textToDoubleList(leafValStr);

    //XXX (debug)
    std::cerr << "\t\tLeaf Values: " << std::endl;
    unsigned int leafIdx = 0;
    
    for (std::vector<double>::iterator iter = leafValues.begin(); iter != leafValues.end(); ++iter)
    {
        std::cerr << "\t\tLeaf " << leafIdx++ << ": " << *iter << std::endl;
    }

    m_failWeight = leafValues[0];
    m_passWeight = leafValues[1];
    std::cerr << "\t\tPass Weight: " << m_passWeight << std::endl;
    std::cerr << "\t\tFail Weight: " << m_failWeight << std::endl;

    uint32_t rectIndex = internalNodes[2];
    tinyxml2::XMLNode *rectNode = rectanglesElement->FirstChild();
    for (uint32_t idx = 0; idx < rectIndex; ++idx)
    {
        rectNode = rectNode->NextSibling();
    }
    m_rectangle = textToUnsignedList(rectNode->FirstChildElement("rect")->GetText());

    std::cerr << "\t\tRectangle: " << std::endl;
    for (std::vector<uint32_t>::iterator iter=m_rectangle.begin(); iter != m_rectangle.end(); ++iter)
    {
        std::cerr <<"\t\t\t" << *iter << " ";
    }
    std::cerr << std::endl;

    // construct lookup table
    for (unsigned int i = 0; i < 8; ++i)
    {
        m_lookupTable[i] = internalNodes[i+3]; 
    }
}

LBPFeature::~LBPFeature()
{
}



// TODO
double LBPFeature::evaluate(IntegralImage &integralImage, 
                            const uint16_t &windowRow, 
                            const uint16_t &windowCol, 
                            const double &scale) const
{
    uint16_t rowStart = windowRow + lround(m_rectangle[1]*scale);
    uint16_t colStart = windowCol + lround(m_rectangle[0]*scale);
    uint16_t rectW = lround(m_rectangle[2]*scale);
    uint16_t rectH = lround(m_rectangle[3]*scale);

    // compute the values for each box
    std::vector<uint32_t> boxVals; 
    for (uint16_t i = 0; i < 3; ++i)
    {
        for(uint16_t j = 0; j < 3; ++j) 
        {
            uint16_t rowMin = rowStart + i*rectH;
            uint16_t rowMax = rowStart + (i+1)*rectH;
            uint16_t colMin = colStart + j*rectW;
            uint16_t colMax = colStart + (j+1)*rectW;
            uint32_t boxUL = integralImage(rowMin,colMin);
            uint32_t boxUR = integralImage(rowMin,colMax);
            uint32_t boxLL = integralImage(rowMax,colMin);
            uint32_t boxLR = integralImage(rowMax,colMax);
            uint32_t boxVal = boxLR + boxUL - boxUR - boxLL;
            boxVals.push_back(boxVal);
        }
    }

    // compare to center box to get local binary pattern
    uint32_t lbp = ((int(boxVals[0]) >= int(boxVals[4])) << 7) | ((int(boxVals[1]) >= int(boxVals[4])) << 6) | ((int(boxVals[2]) >= int(boxVals[4])) << 5) | ((int(boxVals[5]) >= int(boxVals[4])) << 4) | ((int(boxVals[8]) >= int(boxVals[4])) << 3) | ((int(boxVals[7]) >= int(boxVals[4])) << 2) | ((int(boxVals[6]) >= int(boxVals[4])) << 1) | ((int(boxVals[3]) >= int(boxVals[4])) << 0);

    // look up match value in the feature's lookup table
    uint8_t matchValue = (m_lookupTable[lbp >> 5] >> (lbp & 31))&1;

    // set feature score based on match value
    double featureScore;
    if (matchValue == 0)
    {
      featureScore = m_passWeight;
    }
    else
    {
      featureScore = m_failWeight;
    }

    return featureScore;
}



std::vector<uint32_t> LBPFeature::textToUnsignedList(const char *text)
{
    std::vector<uint32_t> unsignedList;
    char *textCpy = strdup(text);
    char **end = &textCpy;
    while (**end != '\0')
    {
        unsignedList.push_back(strtol(textCpy,end,10));
        textCpy = *end;
    }
    return unsignedList;
}


std::vector<double> LBPFeature::textToDoubleList(const char *text)
{
    std::vector<double> doubleList;
    char *textCpy = strdup(text);
    char **end = &textCpy;
    while (**end != '\0')
    {
        doubleList.push_back(strtod(textCpy,end));
        textCpy = *end;
    }
    return doubleList;
}


CascadeStage::CascadeStage(tinyxml2::XMLElement *stageElement, tinyxml2::XMLElement *rectanglesElement)
{
    stageElement->FirstChildElement("stageThreshold")->QueryDoubleText(&m_threshold);

    //XXX (debug)
    std::cerr << "\tThreshold = " << m_threshold << std::endl;
    unsigned int idx = 0;

    // load all features for this stage
    tinyxml2::XMLElement *features = stageElement->FirstChildElement("weakClassifiers");
    for (tinyxml2::XMLElement *featureElement = features->FirstChildElement("_");
         featureElement != NULL; 
         featureElement = featureElement->NextSiblingElement("_"))
    {
        //XXX (debug)
        std::cerr << "\tFeature " << idx++ << std::endl;

        m_features.push_back(LBPFeature(featureElement, rectanglesElement));
    }

    return;
}


CascadeStage::~CascadeStage()
{
}



bool CascadeStage::evaluate(IntegralImage &integralImage, 
                              const uint16_t &row, 
                              const uint16_t &col, 
                              const double &scale) const
{
    double score = 0.0;
    for (std::vector<LBPFeature>::const_iterator feature = m_features.begin(); feature != m_features.end(); ++feature)
    {
        score += feature->evaluate(integralImage, row, col, scale); 
        if (score < m_threshold)
        {
//            std::cerr << "STAGE evaluates to FALSE - score of " << score << " is less than threshold " << m_threshold << std::endl; //XXX
            return false;
        }
    }
//    std::cerr << "STAGE evaluates to TRUE" << std::endl; //XXX
    return true;
}


CascadeClassifier::CascadeClassifier(const std::string &cascadeFile)
{
    // load cascade from xml
    tinyxml2::XMLDocument cascadeDoc;
    tinyxml2::XMLError eResult = cascadeDoc.LoadFile(cascadeFile.c_str());
    tinyxml2::XMLElement *cascade = cascadeDoc.RootElement()->FirstChildElement("cascade");
    
    // parse out general classifier info
    cascade->FirstChildElement("height")->QueryUnsignedText(&m_windowHeight);
    cascade->FirstChildElement("width")->QueryUnsignedText(&m_windowWidth);

    //XXX (debug)
    std::cerr << "Captured from XML: " << std::endl;
    std::cerr << "Window Height = " << m_windowHeight << std::endl;
    std::cerr << "Window Width = " << m_windowWidth << std::endl;
    unsigned int idx = 0;


    // get a pointer to the rectangles
    tinyxml2::XMLElement *rectanglesElement = cascade->FirstChildElement("features");

    // create each of the stages
    tinyxml2::XMLElement *stages = cascade->FirstChildElement("stages");
    for (tinyxml2::XMLElement *stageElement = stages->FirstChildElement("_");
         stageElement != NULL; 
         stageElement = stageElement->NextSiblingElement("_"))
    {
        //XXX (debug)
        std::cerr << "Stage " << idx++ << std::endl;

        m_stages.push_back(CascadeStage(stageElement,rectanglesElement)); 
    }
}


CascadeClassifier::~CascadeClassifier()
{
}


std::vector<detectionLocation> CascadeClassifier::detectMultiScale(const QImage &image, 
                                                                   const double &scaleFactor, 
                                                                   const uint16_t &stepSize) const
{
    // Create integral image from raw QImage
//    std::cerr << "Generating Integral Image" << std::endl;
    IntegralImage integralImage(image);
//    std::cerr << "...DONE!" << std::endl;
    return detectMultiScale(integralImage, scaleFactor, stepSize);
}



std::vector<detectionLocation> CascadeClassifier::detectMultiScale(IntegralImage &integralImage, 
                                                                   const double &scaleFactor, 
                                                                   const uint16_t &stepSize) const
{
    std::vector<detectionLocation> detections;
    double maxScale = std::min(integralImage.width()/m_windowWidth, integralImage.height()/m_windowHeight);
    for (double scale = 1.0; scale < maxScale; scale *= scaleFactor)
    {
        //TODO pass reference to vector instead of returning & appending?
        std::vector<detectionLocation> scaleDetections = detectSingleScale(integralImage, scale, stepSize);
        for (std::vector<detectionLocation>::const_iterator detectionLoc = scaleDetections.begin(); 
             detectionLoc != scaleDetections.end(); ++detectionLoc)
        {
            detections.push_back(*detectionLoc);
        }
    }
    return detections;
}


//TODO carefully review rounding 
std::vector<detectionLocation> CascadeClassifier::detectSingleScale(IntegralImage &integralImage, 
                                                                    const double &scale, 
                                                                    const uint16_t &stepSize) const
{
    std::vector<detectionLocation> detections;
    //XXX subtracting 2 from each dimension of the scan to avoid overflowing past the edge of the integralImage due to rounding
    for (uint16_t windowRow = 0; windowRow < integralImage.height()-lround(m_windowHeight*scale)-2; windowRow += stepSize)
    {
        for (uint16_t windowCol = 0; windowCol < integralImage.width()-lround(m_windowWidth*scale)-2; windowCol += stepSize)
        {
            if (detectAtLocation(integralImage, windowRow, windowCol, scale))
            {
                // detection found
                detectionLocation loc(windowCol, windowRow, lround(m_windowHeight*scale), lround(m_windowWidth*scale));
                detections.push_back(loc);
            }
        }
    }
    return detections;
}


bool CascadeClassifier::detectAtLocation(IntegralImage &integralImage, 
                                         const uint16_t &row, 
                                         const uint16_t &col, 
                                         const double &scale) const
{
    for (std::vector<CascadeStage>::const_iterator stage = m_stages.begin(); stage != m_stages.end(); ++stage)
    {
        if (stage->evaluate(integralImage, row, col, scale) == false) 
        {
            return false;
        }
    }
    return true;
}
