#include "my_belle.h"

int code(int charg, int barion_num, int chanel){
    auto boool = [](int num) { return (num < 0) ? 1 : 0; };
    int code = 0;
    
    code = code + boool(charg) * 1e10;
    code = code + abs(charg) * 1e9;
    
    code = code + boool(barion_num) * 1e8;
    code = code + abs(barion_num) * 1e8;

    code = code + abs(chanel) * 1e4;
    
    return code;
}


void decode(int code, int& charg, int& barion_num, int& chanel) {
    int charg_sign = (code / static_cast<int>(1e10)) % 2;  // Extract sign of charg
    charg = (code / static_cast<int>(1e9)) % static_cast<int>(1e1);  // Extract charg value

    int barion_sign = (code / static_cast<int>(1e8)) % 2;  // Extract sign of barion_num
    barion_num = (code / static_cast<int>(1e7)) % static_cast<int>(1e1);  // Extract barion_num value

    chanel = code % static_cast<int>(1e4);  // Extract chanel value

    // Apply sign to charg and barion_num
    charg = (charg_sign == 1) ? -charg : charg;
    barion_num = (barion_sign == 1) ? -barion_num : barion_num;
}

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1, 
                std::vector<Particle>& child_2, std::vector<Particle>& child_3, std::vector<Particle>& child_4, 
                int charg, int barion_num, int chanel) {
    combination(comb_particle, object, child_1, child_2, child_3, child_4);
    setUserInfo(comb_particle, code(charg, barion_num, chanel));
}

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1,
                int charg, int barion_num, int chanel) {
    combination(comb_particle, object, child_1);
    setUserInfo(comb_particle, code(charg, barion_num, chanel));
}

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1,
                std::vector<Particle>& child_2, int charg, int barion_num, int chanel) {
    combination(comb_particle, object, child_1, child_2);
    setUserInfo(comb_particle, code(charg, barion_num, chanel));
}

template<typename T>
void make_combo(std::vector<Particle>& comb_particle, T object, std::vector<Particle>& child_1,
                std::vector<Particle>& child_2, std::vector<Particle>& child_3,
                int charg, int barion_num, int chanel) {
    combination(comb_particle, object, child_1, child_2, child_3);
    setUserInfo(comb_particle, code(charg, barion_num, chanel));
}
