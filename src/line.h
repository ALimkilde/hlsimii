#pragma once

#include <vector>

struct Webbing {
    double stretch_pct; // Stretch percentage at the tension
    double tension_kN;  // The tension for which the stretch is given in [kN]
    double weight_g_m;  // Weight in [g/m]

    double rho() const; // Linear density [kg/m]
    double ea() const;  // Axial stiffness EA [N]
};

struct Segment {
   Segment(Webbing webbing_main, Webbing webbing_backup, double l_main, double l_backup);

   double ea_main;   // axial stiffness. spring constant times l
   double ea_backup; // axial stiffness. spring constant times l
   double L_main;
   double L_backup;
   double weight;
};

// A highline rig as it is in the bag
class Line {
public:
    explicit Line(std::vector<Segment> segments);

    int num_segments() const;
    const std::vector<Segment>& segments() const;

private:
    std::vector<Segment> segments_;
};
