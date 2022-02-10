#include <iostream>
#include "../include/Settings.h"
#include "../include/Drivingforcesolver.h"
#include "../include/Phasefieldsolver.h"
#include "../include/Solutefieldsolver.h"
#include "../../libtorch/include/torch/script.h"
#include "../include/Temperaturefieldsolver.h"
#include <cmath>

using namespace std;

Drivingforcesolver::Drivingforcesolver(const Settings& localsettings)
{
    this->Load_pytorch_module(localsettings);
    this->Create_fields(localsettings);
}

void Drivingforcesolver::Create_fields(const Settings& localsettings)
{
    this->drivingforce_field.Resize(localsettings._Nx,localsettings._Ny);
}

void Drivingforcesolver::Load_pytorch_module(const Settings& localsettings)
{
    try {
        module = torch::jit::load(localsettings.Pytorch_module_path);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model\n";
        exit(-1);
    }
}

inline double Normalize(double input_val, double mean_val, double std_val)
{
    return (input_val-mean_val)/std_val;
}

void Drivingforcesolver::Calculate_c_phase(const Settings& localsettings, const Phasefieldsolver& local_phasefield, 
                       Solutefieldsolver& localsolutefield, const Temperaturefieldsolver& local_temp)
{
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor input_tensor = torch::ones({1,4});
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        {
            if (local_phasefield.interface_flag[i][j] == 0.5)
            {
                input_tensor[0][0] = Normalize(local_temp.temperature_field[i][j], 876.002902, 32.62102983);//Mean value and mean square deviation of machine learning training set
                input_tensor[0][1] = Normalize(localsolutefield.c_overall_old[i][j], 0.169993173, 0.09522574);
                input_tensor[0][2] = Normalize(1-local_phasefield.phi_old[i][j],0.499914659,0.31618727);
                input_tensor[0][3] = Normalize(local_phasefield.phi_old[i][j],0.500085341,0.31618727);
                inputs.push_back(input_tensor);
                auto output = this->module.forward(inputs).toTensor();

                localsolutefield.c_phase_1[i][j] = output[0][0].item().to<double>();
                localsolutefield.c_phase_2[i][j] = output[0][1].item().to<double>();

                if (localsolutefield.c_phase_1[i][j] <= 0 || localsolutefield.c_phase_1[i][j] >= 1)
                {
                    cout << "c_phase_1 value error! " << localsolutefield.c_phase_1[i][j] << " Location:" << i << " " << j << endl; 
                    cout << local_temp.temperature_field[i][j]<< " "<< localsolutefield.c_overall_old[i][j]<< " "<< local_phasefield.phi_old[i][j] << endl;
                }
                if (localsolutefield.c_phase_2[i][j] <= 0 || localsolutefield.c_phase_2[i][j] >= 1)
                {
                    cout<<"c_phase_2 value error! " << localsolutefield.c_phase_2[i][j] << " Location:" << i << " " << j << endl; 
                }
                inputs.pop_back();
            }
            else if (local_phasefield.interface_flag[i][j] == 0)//liquid 
            {
                localsolutefield.c_phase_1[i][j] = localsolutefield.c_overall_old[i][j];
                localsolutefield.c_phase_2[i][j] = 0;
            }
            else if (local_phasefield.interface_flag[i][j]==1)//solid 
            {
                localsolutefield.c_phase_2[i][j] = localsolutefield.c_overall_old[i][j];
                localsolutefield.c_phase_1[i][j] = 0;
            }
        }
}

inline double G_L(double c_Cu,double T)
{
    if(c_Cu<=0.0 || c_Cu>=1.0)
    {
        cout << "c_Cu value error!" << endl;
        exit(-1);
    } 
    double GHSER_Al = -11276.24+223.048446*T-38.5844296*T*log(T)+18.531982e-3*pow(T,2)-5.764227e-6*pow(T,3)+74092.0/T;
    double GLIQ_Al = 11005.029-11.841867*T+7.934e-20*pow(T,7)+GHSER_Al;
    double GHSER_Cu = -7770.458+130.485235*T-24.112392*T*log(T)-2.65684e-3*pow(T,2)+0.129223e-6*pow(T,3)+52478.0/T;
    double GLIQ_Cu = 12964.736-9.511904*T+5.849e-21*pow(T,7)+GHSER_Cu;
    double G_0 = (1.0-c_Cu)*GLIQ_Al + c_Cu*GLIQ_Cu;
    double L0 = -66622.0+8.1*T;
    double L1 = 46800.0-90.8*T+10.0*T*log(T);
    double L2 = -2812.0;
    double G_id=8.314*T*((1.0-c_Cu)*log(1.0-c_Cu)+c_Cu*log(c_Cu));
    double G_xs=c_Cu*(1.0-c_Cu)*(L0 + L1*(1.0-2.0*c_Cu) + L2*(1.0-2.0*c_Cu)*(1.0-2.0*c_Cu));
    return G_0 + G_id + G_xs;
}

inline double d_G_L(double c_Cu, double T)
{
    if(c_Cu<=0.0 || c_Cu>=1.0)
    {
        cout << "c_Cu value error!" << endl;
        exit(-1);
    } 
    double GHSER_Al = -11276.24+223.048446*T-38.5844296*T*log(T)+18.531982e-3*pow(T,2)-5.764227e-6*pow(T,3)+74092.0/T;
    double GLIQ_Al = 11005.029-11.841867*T+7.934e-20*pow(T,7)+GHSER_Al;
    double GHSER_Cu = -7770.458+130.485235*T-24.112392*T*log(T)-2.65684e-3*pow(T,2)+0.129223e-6*pow(T,3)+52478.0/T;
    double GLIQ_Cu = 12964.736-9.511904*T+5.849e-21*pow(T,7)+GHSER_Cu;
    double d_G_0 = -GLIQ_Al + GLIQ_Cu;
    double L0 = -66622.0+8.1*T;
    double L1 = 46800.0-90.8*T+10.0*T*log(T);
    double L2 = -2812.0;
    double d_G_id=8.314*T*(-log(1.0-c_Cu)+log(c_Cu));
    double d_G_xs=(1.0-2.0*c_Cu)*(L0 + L1*(1.0-2.0*c_Cu) + L2*(1.0-2.0*c_Cu)*(1.0-2.0*c_Cu)) +c_Cu*(1.0-c_Cu)*(-2.0*L1 - L2*4.0*(1.0-2.0*c_Cu));
    return d_G_0 + d_G_id + d_G_xs;
}

inline double G_Al(double c_Cu, double T)
{
    if(c_Cu<=0.0 || c_Cu>=1.0)
    {
        cout << "c_Al value error!" << endl;
        exit(-1);
    }     
    double GHSER_Al = -11276.24+223.048446*T-38.5844296*T*log(T)+18.531982e-3*pow(T,2)-5.764227e-6*pow(T,3)+74092.0/T;
    double GLIQ_Al = 11005.029-11.841867*T+7.934e-20*pow(T,7)+GHSER_Al;
    double GHSER_Cu = -7770.458+130.485235*T-24.112392*T*log(T)-2.65684e-3*pow(T,2)+0.129223e-6*pow(T,3)+52478.0/T;
    double GLIQ_Cu = 12964.736-9.511904*T+5.849e-21*pow(T,7)+GHSER_Cu;
    double G_0 = (1.0-c_Cu)*GHSER_Al + c_Cu*GHSER_Cu;
    double L0 = -53520.0+2.0*T;
    double L1 = 38590.0-2.0*T;
    double L2 = 1170.0;
    double G_id=8.314*T*((1.0-c_Cu)*log(1.0-c_Cu)+c_Cu*log(c_Cu));
    double G_xs=c_Cu*(1.0-c_Cu)*(L0 + L1*(1.0-2.0*c_Cu) + L2*(1.0-2.0*c_Cu)*(1.0-2.0*c_Cu));
    return G_0 + G_id + G_xs;
}

inline double d_G_Al(double c_Cu, double T)
{
    if (c_Cu <= 0.0 || c_Cu >= 1.0)
    {
        cout << "c_Al value error!" << endl;
        exit(-1);
    }
    double GHSER_Al = -11276.24 + 223.048446 * T - 38.5844296 * T * log(T) + 18.531982e-3 * pow(T, 2) - 5.764227e-6 * pow(T, 3) + 74092.0 / T;
    double GLIQ_Al = 11005.029 - 11.841867 * T + 7.934e-20 * pow(T, 7) + GHSER_Al;
    double GHSER_Cu = -7770.458 + 130.485235 * T - 24.112392 * T * log(T) - 2.65684e-3 * pow(T, 2) + 0.129223e-6 * pow(T, 3) + 52478.0 / T;
    double GLIQ_Cu = 12964.736 - 9.511904 * T + 5.849e-21 * pow(T, 7) + GHSER_Cu;
    double d_G_0 = -GHSER_Al + GHSER_Cu;
    double L0 = -53520.0 + 2.0 * T;
    double L1 = 38590.0 - 2.0 * T;
    double L2 = 1170.0;
    double d_G_id = 8.314 * T * (-log(1.0 - c_Cu) + log(c_Cu));
    double d_G_xs = (1.0 - 2.0 * c_Cu) * (L0 + L1 * (1.0 - 2.0 * c_Cu) + L2 * (1.0 - 2.0 * c_Cu) * (1.0 - 2.0 * c_Cu)) + c_Cu * (1.0 - c_Cu) * (-2.0 * L1 - L2 * 4.0 * (1.0 - 2.0 * c_Cu));
    return d_G_0 + d_G_id + d_G_xs;
}

void Drivingforcesolver::Calculate_driving_force(const Settings& localsettings, const Phasefieldsolver& local_phasefield, 
                        const Solutefieldsolver& localsolutefield, const Temperaturefieldsolver& local_temp)
{
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        {
            if (local_phasefield.interface_flag[i][j] == 0.5)
            {
                double c_Cu_L = localsolutefield.c_phase_1[i][j];
                double c_Cu_S = localsolutefield.c_phase_2[i][j];
                double T = local_temp.temperature_field[i][j];
                double d_G = d_G_L(c_Cu_L, T);
                this->drivingforce_field[i][j] = G_L(c_Cu_L, T) - G_Al(c_Cu_S, T) - d_G * (c_Cu_L - c_Cu_S);//Thermodynamic driving force calculation
            }
            else
            {
                this->drivingforce_field[i][j] = 0.0;
            }            
        }
}