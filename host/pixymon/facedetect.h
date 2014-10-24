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

#ifndef FACEDETECT_H 
#define FACEDETECT_H
#include <tinyxml2.h>
#include <stdint.h>
#include <QImage>


class IntegralImage
{
public:
    IntegralImage(const QImage &image);
    ~IntegralImage();

    const uint32_t operator()(const uint16_t &row, const uint16_t &col) const;
    uint32_t operator()(const uint16_t &row, const uint16_t &col);
    uint16_t width();
    uint16_t height();

private:
    uint16_t m_width;
    uint16_t m_height;
    uint32_t m_data[321][201]; //TODO correct height
};


struct detectionLocation
{
    detectionLocation(const uint16_t &x, const uint16_t &y, const uint16_t &w, const uint16_t &h);

    uint16_t locationX;
    uint16_t locationY;
    uint16_t width;
    uint16_t height;
};


// TODO make base feature class and inherit from this
class LBPFeature
{
public:
    LBPFeature(tinyxml2::XMLElement *featureElement, tinyxml2::XMLElement *rectanglesElement);
    ~LBPFeature();

    double evaluate(IntegralImage &integralImage, const uint16_t &windowCol, const uint16_t &windowRow, const double &scale) const;

private:
    std::vector<uint32_t> m_rectangle;
    double m_failWeight;
    double m_passWeight;
    uint32_t m_lookupTable[8];

    std::vector<uint32_t> textToUnsignedList(const char *text);
    std::vector<double> textToDoubleList(const char *text);
};


class CascadeStage
{
public:
    CascadeStage(tinyxml2::XMLElement *stageElement, tinyxml2::XMLElement *rectanglesElement);
    ~CascadeStage();

    bool evaluate(IntegralImage &integralImage, const uint16_t &locationX, const uint16_t &locationY, const double &scale) const;

private:
    double m_threshold;
    std::vector<LBPFeature> m_features; // TODO make this vector of generic feature type
};


class CascadeClassifier
{
public:
    CascadeClassifier(const std::string &xmlFile);
    ~CascadeClassifier();

    std::vector<detectionLocation> detectMultiScale(const QImage& image, const double &scaleFactor, const uint16_t &stepSize) const;

private:
    std::vector<detectionLocation> detectMultiScale(IntegralImage &integralImage, const double &scaleFactor, const uint16_t &stepSize) const;
    std::vector<detectionLocation> detectSingleScale(IntegralImage &integralImage, const double  &scale, const uint16_t &stepSize) const;
    bool detectAtLocation(IntegralImage &integralImage, const uint16_t &locationX, const uint16_t &locationY, const double &scale) const;

    std::vector<CascadeStage> m_stages;
    unsigned int m_windowWidth;
    unsigned int m_windowHeight;
};

#endif // FACEDETECT_H

