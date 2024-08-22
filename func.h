#ifndef MY_FUNC
#define MY_FUNC

#include "my_belle.h"

int code(int charg, int barion_num, int chanel);
int decoder(int code, int& charg, int& barion_num, int& chanel);

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1,
                std::vector<Particle>& child_2, std::vector<Particle>& child_3,
                int charg, int barion_num, int chanel);

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1,
                std::vector<Particle>& child_2, int charg, int barion_num, int chanel);

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1,
                int charg, int barion_num, int chanel);

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1,
                std::vector<Particle>& child_2, std::vector<Particle>& child_3, std::vector<Particle>& child_4,
                int charg, int barion_num, int chanel);

#endif 
