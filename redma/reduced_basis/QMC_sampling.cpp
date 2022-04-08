//
// Created by Federico Betti on 03/04/2022.
//

#include "QMC_sampling.h"

namespace RedMA
{
    QMC_sampling::
    QMC_sampling(unsigned int numSamples)
    {
        M_numSamples = numSamples;

        const double maxAngle = 0.4;
        const double maxAmplitude = 0.2;
        const double maxWidth = 0.3;
        const double minFlow = 0.2;
        const double maxFlow = 0.8;

        // REGISTER ALL POSSIBLE PARAMETERS POTENTIALLY SUBJECT TO SAMPLING

        M_parametersHandler.registerParameter("in1_alphax", 0,
                                              -maxAngle, maxAngle, true, false);
        M_parametersHandler.registerParameter("in1_alphay", 0,
                                              -maxAngle, maxAngle, true, false);
        M_parametersHandler.registerParameter("in1_alphaz", 0,
                                              -maxAngle, maxAngle, true, false);
        M_parametersHandler.registerParameter("in2_alphax", 0,
                                              -maxAngle, maxAngle, true, false);
        M_parametersHandler.registerParameter("in2_alphay", 0,
                                              -maxAngle, maxAngle, true, false);
        M_parametersHandler.registerParameter("in2_alphaz", 0,
                                              -maxAngle, maxAngle, true, false);
        M_parametersHandler.registerParameter("stenosis_amplitude", 0,
                                              0, maxAmplitude, true, false);
        M_parametersHandler.registerParameter("stenosis_width", 0,
                                              0, maxWidth, true, false);
        M_parametersHandler.registerParameter("flow_rate", 0,
                                              minFlow, maxFlow, true, false);

        std::tie(M_paramsNames, M_paramsBounds) = setParametersToBeSampled(M_parametersHandler);

        // while there is no check here, to achieve the optimal bound on low discrepancy one should pick
        // all the entries of the generating vector such that they don't have factors in common with M_numSamples
        M_generatingVector = {3, 7, 11, 13, 17, 19, 23, 29, 31};

        // checking consistency
        checkConsistency(M_paramsBounds, M_generatingVector);

        // M_samplesVector = getSamples(M_numSamples, M_generatingVector, M_paramsVector, M_paramsBounds);
    }

    std::tuple<std::vector<std::string>, std::vector<std::vector<double>>>
    QMC_sampling::
    setParametersToBeSampled(GeometricParametersHandler& paramsHandler)
    {
        std::vector<std::string> paramsNames;
        std::vector<std::vector<double>> paramsBounds;
        for (auto it = paramsHandler.getParametersMap().begin();
                            it != paramsHandler.getParametersMap().end(); ++ it)
        {
                 paramsNames.push_back(it->first);
                 paramsBounds.push_back(GetBounds(it->first, paramsHandler));
        }
        return std::make_tuple(paramsNames, paramsBounds);
    }

    void
    QMC_sampling::
    checkConsistency(std::vector<std::vector<double>> paramsBounds, std::vector<double> generatingVector)
    {
        if (paramsBounds.size() != generatingVector.size())
        {
            throw new Exception("The dimension of the parameters vector doesn't correspond with"
                            "the one of the generating vector");
        }
    }

    std::vector<double>
    QMC_sampling::
    GetBounds(std::string paramName, GeometricParametersHandler paramsHandler)
    {
        std::vector<double> bounds;
        if (paramsHandler.exists(paramName))
        {
            bounds.push_back(paramsHandler.getParametersMap().at(paramName)->getMinValue());
            bounds.push_back(paramsHandler.getParametersMap().at(paramName)->getMaxValue());
        }
        return bounds;
    }

    std::vector<double>
    QMC_sampling::
    MultiplyVectorByScalar(std::vector<double> v, double k)
    {
    std::vector<double> v1(v.size());
    std::transform(v.begin(), v.end(), v1.begin(),
                   [k](auto& c) {return c * k;});
    return v1;
    }

    std::vector<double>
    QMC_sampling::
    getFractionalPart(std::vector<double> v)
    {
        std::transform(v.begin(), v.end(), v.begin(),
                [](double &c){return c-std::floor(c); });
        return v;
    }

    std::vector<double>
    QMC_sampling::
    translateSampleToBounds(std::vector<double> sample, std::vector<std::vector<double>> paramsBounds)
    {
        for (unsigned int i = 0; i < sample.size(); i++)
        {
            sample[i] = paramsBounds[i][0] + sample[i] * (paramsBounds[i][1] - paramsBounds[i][0]);
        }
        return sample;
    }

    std::vector<std::map<std::string, double>>
    QMC_sampling::
    getSamples(unsigned int N, std::vector<double> generatingVector,
               std::vector<std::string> paramsNames,  std::vector<std::vector<double>> paramsBounds)
    {
        std::vector<std::map<std::string, double>> samplesVector;
        unsigned int d = paramsNames.size();
        for (unsigned int i = 0; i < N; i ++)
        {
            std::vector<double> current_sample = MultiplyVectorByScalar(generatingVector, i/double(N));
            current_sample = getFractionalPart(current_sample);
            current_sample = translateSampleToBounds(current_sample, paramsBounds);
            std::map<std::string, double> map;
            for (unsigned int j = 0; j < d; ++ j)
            {
                map.insert(std::pair<std::string, double> (paramsNames[j], current_sample[j]));
            }
            samplesVector.push_back(map);
        }
        return samplesVector;
    }

}

