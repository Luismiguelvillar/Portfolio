/*
Final Project c++
-----------------
Simulation of particle detection using
Atlas-like detection system. Code runs strictly in keV units for electrons.

by Luis Miguel Villar Padruno 
30/03/2024
*/
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <random>

using std::cout;
using std::endl;
using std::cin;
using std::string;
using std::vector;

// The efficiencies run from 0 to 10
// Each layer of the tracker and muon chamber have a chance of detecting given by these values.
// For the calorimeter calorimeter_efficiency represents the chance of correct functioning (detection).
// The true efficiencies (the action of multilayered trackers / muons) will be calculated using this starting values.
const double detector_analysis_efficiency{10};
const double tracker_efficiency{10};
const double muon_chamber_efficiency{10};
const double likelihood_leptonic_decay{5};
const double calorimeter_efficiency{10};
int ids{0};
int seed;

// Comentarios 
// 1. Se debe mejorar la print function de todas las clases.
// 2. Probar el contador de energia con varios eventos diferentes. Luego pensar en ya terminar.

// Mathematical functions
int number_of_detections_for_confidence(double detector_efficiency, int number_of_detectors)
{
  // Here we will use the bonomial distribution to model how confident we are in a result. And more importantely, how many
  // detector of a total number of detectors must output TRUE to consider a detection reasonable for a 95% confidence
  // level. Here the efficiency and the desired confidence level are number from 0 to 1.
  double result{0};
  double expected_value{number_of_detectors * detector_efficiency};
  double standard_deviation{sqrt(expected_value*(1 - detector_efficiency))};

  // We will be 95% secure of our result when the number of detector detecting the outcome is within 2 sigmas of the mean.
  result = expected_value - 2*standard_deviation;
  return result;
}

vector<double> divide_random(double energy, int parts)
{
  // This fucntions takes the energy and divides it in random parts that add up to the initial energy.
  double original_energy{energy};
  static double errors{0}; // using static variable
  vector<double> divisions(parts);
  for(int i{0}; i < parts - 1; i++)
  {
    if((int)energy - parts + 1 > 0) 
    {
     divisions[i] = (float)(1 + rand() % ((int)energy - parts + 1));
    }
    else // If we get a division by zero we just devide equally an avoid the error (really rare error, around 0.3 % of the times for Neutrons, for example)
    {
      errors = errors + 1;
      double energy_per_division = original_energy / (double)parts;
      for(int i{0}; i < parts; i++)
      {
        divisions[i] = energy_per_division;
      }
      return divisions;
    }
    energy = energy - divisions[i];
  }
  divisions.back() = energy;
  return divisions;
}

// ---Classes---
//FOUR MOMENTUM
class FourMomentum 
{
  protected:
    double energy_over_c{0};
    double x1_coordinate{0};
    double x2_coordinate{0};
    double x3_coordinate{0};
  public:
    FourMomentum() = default;
    FourMomentum(double energy_over_c_initial, double x1_initial, double x2_initial, double x3_initial)
    {
      if(energy_over_c_initial < 0)
      {
        cout << "Invalid value for the Energy/c: ";
        cin>>energy_over_c_initial;
        while (cin.fail() or energy_over_c_initial < 0)
        {
          cout<<"Invalid Energy/c, try again: ";
          cin>>energy_over_c_initial;
        }
        energy_over_c = energy_over_c_initial;
      } else
      {
        energy_over_c = energy_over_c_initial;
      }
      x1_coordinate = x1_initial;
      x2_coordinate = x2_initial;
      x3_coordinate = x3_initial;
    }
    ~FourMomentum(){cout<<"Four Momentum vector destroyed"<<endl;}
    void print_FourMomentum();
    // Copy Constructor
    FourMomentum(FourMomentum const & object)
    {
      energy_over_c = object.energy_over_c;
      x1_coordinate = object.x1_coordinate;
      x2_coordinate = object.x2_coordinate;
      x3_coordinate = object.x3_coordinate;
    }
    //Move Constructor
    FourMomentum( FourMomentum && object)
    {
      energy_over_c = object.energy_over_c;
      x1_coordinate = object.x1_coordinate;
      x2_coordinate = object.x2_coordinate;
      x3_coordinate = object.x3_coordinate;
      //Now we clear the memory of the last vector
      energy_over_c = 0;
      x1_coordinate = 0;
      x2_coordinate = 0;
      x3_coordinate = 0;
    }
    //Setters
    void set_energy_over_c(const double energy_over_c_set) 
    {
      if(energy_over_c_set >= 0)
      {
        energy_over_c =  energy_over_c_set;
      } else
      {
        cout<<"Energy is negative, the absolute value will be set."<<endl;
        energy_over_c = (-1)*energy_over_c_set;
      }
    }
    void set_x1(const double x1_set) {x1_coordinate = x1_set;}
    void set_x2(const double x2_set) {x2_coordinate = x2_set;}
    void set_x3(const double x3_set) {x3_coordinate = x3_set;}
    //Getters
    double get_energy_over_c() const {return energy_over_c;}
    double get_x1() const {return x1_coordinate;}
    double get_x2() const {return x2_coordinate;}
    double get_x3() const {return x3_coordinate;}

    //Overloading Operators & Dot Product Function
    FourMomentum operator+(FourMomentum const & vector_sumado) // The action of this is to change the coordinates of the first vector in the operation v1 + v2
    {
      energy_over_c = energy_over_c + vector_sumado.get_energy_over_c();
      x1_coordinate = x1_coordinate + vector_sumado.get_x1();
      x2_coordinate = x2_coordinate + vector_sumado.get_x2();
      x3_coordinate = x3_coordinate + vector_sumado.get_x3();
      return *this;
    }
    double dot_product(FourMomentum const & vector_dot)
    {
      return (energy_over_c*vector_dot.get_energy_over_c() - (x1_coordinate*vector_dot.get_x1() + x2_coordinate*vector_dot.get_x2() + x3_coordinate*vector_dot.get_x3()));
    }
    //Assignment operator
    FourMomentum & operator=(FourMomentum const & object)
    {
      if(&object == this) return *this;
      energy_over_c = object.energy_over_c;
      x1_coordinate = object.x1_coordinate;
      x2_coordinate = object.x2_coordinate;
      x3_coordinate = object.x3_coordinate;
      return *this;
    }
    //Move Assignment operator
    FourMomentum & operator=(FourMomentum && object_to_migrate)
    {
      std::swap(energy_over_c, object_to_migrate.energy_over_c);
      std::swap(x1_coordinate, object_to_migrate.x1_coordinate);
      std::swap(x2_coordinate, object_to_migrate.x2_coordinate);
      std::swap(x3_coordinate, object_to_migrate.x3_coordinate);
      return *this;
    }

}; // End of Four Momentum Class
void FourMomentum::print_FourMomentum()
{
  cout<<"Four Momentum = ("<<energy_over_c<<", "<<x1_coordinate<<", "<<x2_coordinate<<", "<<x3_coordinate<<")"<<endl;
}
//PARTICLE
class Particle
{
  protected:
    unsigned int identifier{0};
    double true_energy{0};
    double spin{0};
    bool is_anti{false};
    vector<int> tracker_seeds{rand() % 10 + 1, rand() % 10 + 1, rand() % 10 + 1};
    vector <int> muon_seed{rand() % 10 + 1, rand() % 10 + 1};
    int calorimeter_seed{rand() % 10 + 1};

  public:
    Particle() 
    {
      identifier = ids; //AQUI
      true_energy = 0;
      spin = 0;
      is_anti = false;
      ids = ids+1; //AQUI
    };
    Particle(int identifier_initial, double true_energy_initial, double spin_initial, bool is_anti_initial)
    {
      identifier = ids; //AQUI
      ids = ids + 1; //AQUI  
      true_energy = true_energy_initial;
      spin = spin_initial;
      is_anti = is_anti_initial;
    }
    virtual ~Particle() {/*cout<<"Particle Destroyed"<<endl;*/}
    // Copy Constructor
    Particle(Particle const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      // Seeds will be kept random.
    }
    // Move Constructor
    Particle(Particle && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      // Now we clear memory
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
    }
  
    // Assignment operator
    Particle & operator=(Particle const & object)
    {
      if(&object == this) return *this;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      // Seeds will be kept random.
      return *this;
    }

    // Move Assignment operator
    Particle & operator=(Particle && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      return *this;
    }
    // Setters 
    virtual void set_identifier(const int id_set) {identifier = id_set;}
    virtual void set_true_energy(const double true_energy_set)
    {
      if(true_energy_set >= 0 )
      {
        true_energy = true_energy_set;
      } else
      {
        cout<<"True energy must be positive. Absolute value will be considered."<<endl;
        true_energy = true_energy_set*(-1);
      }
    } 
    virtual void set_spin(const double spin_set)
    {
      if(spin > 0)
      {
        spin = spin_set;
      } else
      {
        cout<<"Spin must be a positive quantity"<<endl;
        spin = spin_set*(-1);
      }
    }
    virtual void set_is_anti(const bool anti_set)
    {
      if(anti_set == false)
      {
        is_anti = false;
      } else
      {
        is_anti = true;
      }
    }
    virtual void change_seeds_tracker()
    {
      tracker_seeds[0] = rand() % 10 + 1;
      tracker_seeds[1] = rand() % 10 + 1;
      tracker_seeds[2] = rand() % 10 + 1; 
    }
    virtual void change_seed_calorimeter() {calorimeter_seed = rand() % 10 + 1;}
    virtual void change_seeds_muon_chamber()
    {
      muon_seed[0] = rand() % 10 + 1;
      muon_seed[1] = rand() % 10 + 1;
    }
  
    // Getters
    virtual int get_identifier() const {return identifier;}
    virtual double get_true_energy() const {return true_energy;}
    virtual double get_spin() const {return spin;}
    virtual int get_calorimeter_seed() const {return calorimeter_seed;}
    virtual bool get_is_anti() const {return is_anti;}
    virtual int get_charge() const {return 0;}
    virtual double get_mass() const {return 0;}
    virtual vector<int> get_tracker_seed() {return tracker_seeds;}
    virtual double energy_layer_n(const int n)
    {
      cout<<"You have reached the Particle Class, something wrong happened..."<<endl;
      return 0;
    }
    virtual vector<int> get_muon_seeds() const {return muon_seed;}
    virtual string get_decay_type() const
    {
      cout<<"Particle does not decay."<<endl;
      return "None";
    }  
    virtual vector<std::shared_ptr<Particle>> get_decay_products()
    {
      cout<<"Particle does not decay"<<endl;
      vector<std::shared_ptr<Particle>> result{std::shared_ptr<Particle>(this)};
      return result; // returns a vector containing a pointer to itself.
    }
    virtual bool get_decayed_yet() const
    {
      cout<<"This particle does not decay. The return will be false"<<endl;
      return false;
    }

    // Decay virtual function
    virtual void decay()
    {
      cout<<"Particle does not decay. Decay function will not run"<<endl;
    }
    // Print Particle
    virtual void print_name()
    {
      cout<<"Particle: [ID, spin, is anti] = ["<<identifier<<", "<<spin<<", "<<is_anti<<"]. ";
      cout<<"The true energy is: "<<true_energy<<endl;
    }
  };

class Lepton : public Particle
{
  protected:
    int charge{0};
    double mass{0};

  public:
    Lepton() = default;
    Lepton(int charge_initial, double mass_initial)
    {
      if(charge_initial == 1 or charge_initial == - 1 or charge_initial== 0)
      {
        charge = charge_initial;
        return;
      }
      bool valid_input{false};
      int new_charge{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid charge for a lepton. Try again: ";
          cin>>new_charge;
          if(new_charge != 1 and new_charge != 0 and new_charge != -1)
          {
            throw std::invalid_argument("----> Charge of lepton must the either 1, 0 or -1 <----");
          } else
          {
            valid_input = true;
            charge = charge_initial;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
    
      if(mass_initial >= 0)
      {
        mass = mass_initial;
        return;
      }
      valid_input = false;
      int new_mass{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid mass for a lepton. Try again: ";
          cin>>new_mass;
          if(new_mass < 0)
          {
            throw std::invalid_argument("----> Mass of lepton must be greater than 0 <----");
          } else
          {
            valid_input = true;
            mass = mass_initial;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
    }
    virtual ~Lepton() {/*cout<<"Lepton Destroyed"<<endl;*/}
    // Copy Constructor
    Lepton(Lepton const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
    }
    // Move Constructor
    Lepton(Lepton && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      // Now we clear
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
    } 
    // Assignment Operator
    Lepton & operator=(Lepton const & object)
    {
      if(&object == this) return *this;
      mass = object.mass;
      charge = object.charge;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      return *this;
    }
    // Move Assignment Operator
    Lepton & operator=(Lepton && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      return *this;
    }

    // Setters
    virtual void set_charge(const int charge_set)
    {
      if(charge_set == 1 or charge_set == - 1 or charge_set == 0)
      {
        charge = charge_set;
        return;
      }
      bool valid_input{false};
      int new_charge{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid charge for a lepton. Try again: ";
          cin>>new_charge;
          if(new_charge != 1 and new_charge != 0 and new_charge != -1)
          {
            throw std::invalid_argument("----> Charge of lepton mus the either 1, 0 or -1 <----");
          } else
          {
            valid_input = true;
            charge = charge_set;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
    }
    virtual void set_mass(const double mass_set)
    {
      if(mass_set >= 0)
      {
        mass = mass_set;
        return;
      }
      bool valid_input{false};
      int new_mass{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid mass for a lepton. Try again: ";
          cin>>new_mass;
          if(new_mass < 0)
          {
            throw std::invalid_argument("----> Mass of lepton must be greater than 0 <----");
          } else
          {
            valid_input = true;
            mass = mass_set;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
    }

    // Getters
    virtual double get_mass() const {return mass;}
    int get_charge() const {return charge;}

    // Print Lepton
    virtual void print_name()
    {
      cout<<"Lepton: [ID, mass, charge, spin, is anti] = "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<"]. ";
      cout<<"The true energy is: "<<true_energy<<endl;
    }

};

class Quark : public Particle 
{
  protected:
  char color_charge{'r'}; // Char can be 'r' & 'g' & 'b'.
  double mass{0};
  double charge{0};

  public:
  Quark() = default;
  Quark(char color_charge_initial, double mass_initial, double charge_initial)
  {
    if(color_charge_initial != 'r' and color_charge_initial != 'g'  and color_charge_initial != 'b')
    {
      cout<<"Color charge must be either 'r', 'g' or 'b':";
      cin>>color_charge_initial;
      while(cin.fail() or (color_charge_initial != 'r' and color_charge_initial != 'g' and color_charge_initial != 'b'))
      {
        cout<<"Invalid input, try again: ";
        cin.clear();
        cin.ignore();
        cin>>color_charge_initial;
      }
      color_charge = color_charge_initial;
    } else
    {
      color_charge = color_charge_initial;
    }

    if(mass_initial >= 0)
    {
      mass = mass_initial;
      return;
    }
    bool valid_input{false};
    int new_mass{0};
    while(valid_input != true)
    {
      try
      {
        cout<<"Invalid mass for a quark. Try again: ";
        cin>>new_mass;
        if(new_mass < 0)
        {
          throw std::invalid_argument("----> Mass of quark must be greater than 0 <----");
        } else
        {
          valid_input = true;
          mass = mass_initial;
        }
      }
      catch(const std::exception& e)
      {
        std::cerr<<"Error: "<<e.what()<<endl;
        cin.clear();
        cin.ignore();
      }
    }

    if(charge_initial == 0.3333 or charge_initial == - 0.3333 or charge_initial== 0.6666 or charge_initial == -0.6666)
    {
      charge = charge_initial;
      return;
    }
    valid_input = false;
    int new_charge{0};
    while(valid_input != true)
    {
      try
      {
        cout<<"Invalid charge for a quark. Try again: ";
        cin>>new_charge;
        if(new_charge != 0.3333 and new_charge != -0.3333 and new_charge != -0.6666 and new_charge != 0.6666)
        {
          throw std::invalid_argument("----> Charge of quark must the either +-1/3 [0.3333], or +-2/3 [0.6666] <----");
        } else
        {
          valid_input = true;
          charge = charge_initial;
        }
      }
      catch(const std::exception& e)
      {
        std::cerr<<"Error: "<<e.what()<<endl;
        cin.clear();
        cin.ignore();
      }
    }
  }
  virtual ~Quark() {cout<<"Quark Destroyed"<<endl;}
  // Copy Constructor 
  Quark(Quark const & object)
  {
    identifier = object.identifier;
    true_energy = object.true_energy;
    spin = object.spin;
    is_anti = object.is_anti;
    color_charge = object.color_charge;
    mass = object.mass;
    charge = object.mass;
  }
  // Move Constructor
  Quark(Quark && object)
  {
    identifier = object.identifier;
    true_energy = object.true_energy;
    spin = object.spin;
    is_anti = object.is_anti;
    charge = object.charge;
    mass = object.mass;
    color_charge = object.color_charge;
    // Now we clear
    charge = 0;
    mass = 0; 
    identifier = -1; // Unconsidered
    true_energy = 0;
    spin = 0;
    is_anti = false;
    tracker_seeds.clear();
    muon_seed.clear();
    calorimeter_seed = 0;
    color_charge = 'n'; // n = none
  }
  // Assignment Operator
  Quark & operator=(Quark const & object)
  { 
    if(&object == this) return *this;
    mass = object.mass;
    charge = object.charge;
    identifier = object.identifier;
    true_energy = object.true_energy;
    spin = object.spin;
    is_anti = object.is_anti;
    color_charge = object.color_charge;
    return *this;
  }
  // Move Assignment Operator
  Quark & operator=(Quark && object)
  {
    std::swap(identifier, object.identifier);
    std::swap(true_energy, object.true_energy);
    std::swap(spin, object.spin);
    std::swap(is_anti, object.is_anti);
    std::swap(mass, object.mass);
    std::swap(charge, object.charge);
    std::swap(color_charge, object.color_charge);
    return *this;
  }

  // Setters
  void set_color_charge(const char color_charge_set)
  {
    if(color_charge_set == 'r' or color_charge_set == 'g' or color_charge_set == 'b')
    {
      color_charge = color_charge_set;
    } else
    {
      cout<<"The value introduced does not represent a valid color charge."<<endl;
      color_charge = color_charge;
    }
  }
  virtual void set_mass(const double mass_set) 
  {
    if(mass_set >= 0)
    {
      mass = mass_set;
      return;
    }
    bool valid_input{false};
    int new_mass{0};
    while(valid_input != true)
    {
      try
      {
        cout<<"Invalid mass for a quark. Try again: ";
        cin>>new_mass;
        if(new_mass < 0)
        {
          throw std::invalid_argument("----> Mass of quark must be greater than 0 <----");
        } else
        {
          valid_input = true;
          mass = mass_set;
        }
      }
      catch(const std::exception& e)
      {
        std::cerr<<"Error: "<<e.what()<<endl;
        cin.clear();
        cin.ignore();
      }
    }
  }
  virtual void set_charge(const double charge_set)
  {
    if(charge_set == 0.3333 or charge_set == - 0.3333 or charge_set == 0.6666 or charge_set == -0.6666)
    {
      charge = charge_set;
      return;
    }
    bool valid_input{false};
    int new_charge{0};
    while(valid_input != true)
    {
      try
      {
        cout<<"Invalid charge for a quark. Try again: ";
        cin>>new_charge;
        if(new_charge != 0.3333 and new_charge != -0.3333 and new_charge != -0.6666 and new_charge != 0.6666)
        {
          throw std::invalid_argument("----> Charge of quark must the either +-1/3 [0.3333], or +-2/3 [0.6666] <----");
        } else
        {
          valid_input = true;
          charge = charge_set;
        }
      }
      catch(const std::exception& e)
      {
        std::cerr<<"Error: "<<e.what()<<endl;
        cin.clear();
        cin.ignore();
      }
    }
  }
  
  // Getters
  virtual double get_mass() {return mass;}
  virtual char get_color_charge() {return color_charge;}
  virtual int get_charge() {return charge;}

  // Print Arbitrary Quark
  virtual void print_name()
  {
    cout<<"Quark: [ID, mass, charge, spin, is anti, color charge] = "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<", "<<color_charge<<"]. ";
    cout<<"The true energy is: "<<true_energy<<endl;
  }

};

class Boson : public Particle 
{
  protected:
    double mass{0};
    int charge{0};

  public:
    Boson() = default;
    Boson(double mass_initial, int charge_initial)
    {
      if(mass_initial >= 0)
      {
        mass = mass_initial;
        return;
      }
      bool valid_input{false};
      int new_mass{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid mass for a boson. Try again: ";
          cin>>new_mass;
          if(new_mass < 0)
          {
            throw std::invalid_argument("----> Mass of boson must be greater than 0 <----");
          } else
          {
            valid_input = true;
            mass = mass_initial;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
      
      if(charge_initial == 1 or charge_initial == - 1 or charge_initial== 0)
      {
        charge = charge_initial;
        return;
      }
      valid_input = false;
      int new_charge{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid charge for a boson. Try again: ";
          cin>>new_charge;
          if(new_charge != 1 and new_charge != 0 and new_charge != -1)
          {
            throw std::invalid_argument("----> Charge of boson must the either 1, 0 or -1 <----");
          } else
          {
            valid_input = true;
            charge = charge_initial;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
    }
    virtual ~Boson() {cout<<"Boson Destroyed"<<endl;}
    Boson(Boson const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
    }
    // Move Constructor
    Boson(Boson && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      // Now we clear
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
    }
    // Assignment Operator
    Boson & operator=(Boson const & object)
    {
      if(&object == this) return *this;
      mass = object.mass;
      charge = object.charge;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      return *this;
    }
    // Move Assignment Operator
    Boson & operator=(Boson && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      return *this;
    }
    
    // Setters
    virtual void set_mass(const double mass_set)
    {
      if(mass_set >= 0)
      {
        mass = mass_set;
        return;
      }
      bool valid_input{false};
      int new_mass{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid mass for a boson. Try again: ";
          cin>>new_mass;
          if(new_mass < 0)
          {
            throw std::invalid_argument("----> Mass of boson must be greater than 0 <----");
          } else
          {
            valid_input = true;
            mass = mass_set;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
    }
    virtual void set_charge(const int charge_set)
    {
      if(charge_set == 1 or charge_set == - 1 or charge_set == 0)
      {
        charge = charge_set;
        return;
      }
      bool valid_input{false};
      int new_charge{0};
      while(valid_input != true)
      {
        try
        {
          cout<<"Invalid charge for a lepton. Try again: ";
          cin>>new_charge;
          if(new_charge != 1 and new_charge != 0 and new_charge != -1)
          {
            throw std::invalid_argument("----> Charge of lepton mus the either 1, 0 or -1 <----");
          } else
          {
            valid_input = true;
            charge = charge_set;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        }
      }
    }
    // Getters
    virtual double get_mass() const {return mass;}
    virtual int get_charge() const {return charge;}
};

class Proton : public Particle
{
  protected:
    int charge{1};
    double mass{938}; // Mev/c^2
    vector<double> energies_layers{0,0,0,0,0,0,0,0,0};

  public:
    Proton()
    {
      charge = 1; mass = 938; true_energy = mass;
      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
      spin = 0.5;
    }
    Proton(double true_energy_initial)
    {
      charge = 1;
      mass = 938;
      spin = 0.5;
      true_energy = mass;
      // Try and catch implementation...
      if(true_energy_initial >= mass)
      {
        true_energy = true_energy_initial;
        vector<double> divisions_true_energy{divide_random(true_energy, 7)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = divisions_true_energy[2]; //HAD1
        energies_layers[3] = divisions_true_energy[3]; //HAD2
        energies_layers[4] = divisions_true_energy[4]; //TRACK1
        energies_layers[5] = divisions_true_energy[5]; //TRACK2
        energies_layers[6] = divisions_true_energy[6]; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2

        return;
      }
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }

      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    ~Proton() {/*cout<<"Proton Destroyed"<<endl;*/}
    // Copy Constructor
    Proton(Proton const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
    }
    // Move Constructor
    Proton(Proton && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      } 
      // Now we clear
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
      energies_layers.clear(); // OJO AQUI
    }
    // Assignment Operator
    Proton & operator=(Proton const & object)
    {
      if(&object == this) return *this;
      mass = object.mass;
      charge = object.charge;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
      return *this;
    }
    // Move Assignment Operator
    Proton & operator=(Proton && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      for(int i{0}; i < 9; i++)
      {
        std::swap(energies_layers[i], object.energies_layers[i]);
      }
      return *this;
    }

    // Setters
    void set_true_energy(const double true_energy_set)
    {
      if(true_energy_set >= mass)
      {
        true_energy = true_energy_set;
        vector<double> divisions_true_energy{divide_random(true_energy, 7)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = divisions_true_energy[2]; //HAD1
        energies_layers[3] = divisions_true_energy[3]; //HAD2
        energies_layers[4] = divisions_true_energy[4]; //TRACK1
        energies_layers[5] = divisions_true_energy[5]; //TRACK2
        energies_layers[6] = divisions_true_energy[6]; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2
        return;
      }
      // Try and catch implementation...
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }
      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }

    // Getters
    double energy_layer_n(const int n)
      {
        if(n == 0 or n == 1 or n == 2 or n == 3 or n == 4 or n == 5 or n== 6 or n == 7 or n == 8)
        {
          return energies_layers[n];
        } else
        {
          cout<<"Index out of bounds. Returning the content of the first layer."<<endl;
          return energies_layers[0];
        }
      }
    int get_charge() const {return charge;}
    double get_mass() const {return mass;}
    // Print Proton
    void print_name()
    {
      cout<<"Proton: [ID, mass, charge, spin, is anti] "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<"]"<<endl;
      cout<<"Expected energy Calorimeter: [EM_1: "<<energies_layers[0]<<", EM_2: "<<energies_layers[1]<<", HAD_1: "<<energies_layers[2]<<", HAD_2: "<<energies_layers[3]<<"]."<<endl;
      cout<<"Expected energy Tracker: [Inner: "<<energies_layers[4]<<", Outer: "<<energies_layers[5]<<", Strip: "<<energies_layers[6]<<"]"<<endl;
      cout<<"Expected energy Muon Chamber: [Inner: "<<energies_layers[7]<<", Outer: "<<energies_layers[8]<<"]"<<endl;
      cout<<"The true energy is: "<<true_energy<<endl;
    }
};

class Neutron : public Particle
{
  protected:
  int charge{0};
  double mass{939.5}; 
  vector<double> energies_layers{0,0,0,0,0,0,0,0,0};

  public:
    Neutron() 
    {
      mass = 939.5;
      charge = 0;
      spin = 0.5;
      true_energy = mass;

      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    Neutron(double true_energy_initial)
    {
      charge = 0;
      spin = 0.5;
      mass = 939.5;
      true_energy = mass;

      // Try and catch implementation...
      if(true_energy_initial >= mass)
      {
        true_energy = true_energy_initial;
        vector<double> divisions_true_energy{divide_random(true_energy, 7)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = divisions_true_energy[2]; //HAD1
        energies_layers[3] = divisions_true_energy[3]; //HAD2
        energies_layers[4] = divisions_true_energy[4]; //TRACK1
        energies_layers[5] = divisions_true_energy[5]; //TRACK2
        energies_layers[6] = divisions_true_energy[6]; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2
        return;
      }
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }

      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
      is_anti =false;
    }
    ~Neutron() {/*cout<<"Neutron Destroyed"<<endl;*/}
    // Copy Constructor
    Neutron(Neutron const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
    }
    // Move Constructor
    Neutron(Neutron && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      } 
      // Now we clear
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
      energies_layers.clear(); // OJO AQUI
    }
    // Assignment Operator
    Neutron & operator=(Neutron const & object)
    {
      if(&object == this) return *this;
      mass = object.mass;
      charge = object.charge;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
      return *this;
    }
    // Move Assignment Operator
    Neutron & operator=(Neutron && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      for(int i{0}; i < 9; i++)
      {
        std::swap(energies_layers[i], object.energies_layers[i]);
      }
      return *this;
    }

    // Setters
    void set_energies_in_layers(const double set_energy0, const double set_energy1, const double set_energy2, const double set_energy3)
      {
        if(set_energy0 + set_energy1 + set_energy2 + set_energy3 == true_energy and (set_energy0 > 0) and (set_energy1 > 0) and (set_energy2 > 0) and (set_energy3 > 0))
        {
          energies_layers[0] = set_energy0;
          energies_layers[1] = set_energy1;
          energies_layers[2] = set_energy2;
          energies_layers[3] = set_energy3;
        } else
        {
          cout<<"The values given does not sum up to the energy or one of them is negative. The energy/c of the electron"<<endl;
          cout<<"will be divided equally between the four layers."<<endl;
          energies_layers[0] = true_energy/4;
          energies_layers[1] = true_energy/4;
          energies_layers[2] = true_energy/4;
          energies_layers[3] = true_energy/4;
        }
      }
    void set_mass() 
    {
      cout<<"The mass of the neutron is 939.5 Mev/c^2"<<endl;
      mass = 939.5*pow(10,6);
    }
    void set_charge()
    {
      cout<<"Neutrons have charge of 0"<<endl;
      charge = 0;
    }
    void set_true_energy(const double true_energy_set)
    {
      if(true_energy_set >= mass)
      {
        true_energy = true_energy_set;
        vector<double> divisions_true_energy{divide_random(true_energy, 7)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = divisions_true_energy[2]; //HAD1
        energies_layers[3] = divisions_true_energy[3]; //HAD2
        energies_layers[4] = divisions_true_energy[4]; //TRACK1
        energies_layers[5] = divisions_true_energy[5]; //TRACK2
        energies_layers[6] = divisions_true_energy[6]; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2
        return;
      }
    // Try and catch implementation...
    double new_energy;
    bool valid_input{false};
    while (valid_input != true)
    {
      try
      {
        cout<<"Try introducing another energy:";
        cin >> new_energy;
        if(new_energy <= mass)
        {
          throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
        } else
        {
          valid_input = true;
          true_energy = new_energy;
        }
      }
      catch(const std::exception& e)
      {
        std::cerr<<"Error: "<<e.what()<<endl;
        cin.clear();
        cin.ignore();
      } 
    }
    vector<double> divisions_true_energy{divide_random(true_energy, 7)};
    energies_layers[0] = divisions_true_energy[0]; //EM1
    energies_layers[1] = divisions_true_energy[1]; //EM2
    energies_layers[2] = divisions_true_energy[2]; //HAD1
    energies_layers[3] = divisions_true_energy[3]; //HAD2
    energies_layers[4] = divisions_true_energy[4]; //TRACK1
    energies_layers[5] = divisions_true_energy[5]; //TRACK2
    energies_layers[6] = divisions_true_energy[6]; //TRACK3
    energies_layers[7] = 0; //MUON1
    energies_layers[8] = 0; //MUON2
  }
    // Getters
    double energy_layer_n(const int n)
      {
        if(n == 0 or n == 1 or n == 2 or n == 3 or n == 4 or n == 5 or n == 6 or n == 7 or n == 8)
        {
          return energies_layers[n];
        } else
        {
          cout<<"Index out of bounds. Returning the content of the first layer."<<endl;
          return energies_layers[0];
        }
      }
    int get_charge() const {return charge;}
    double get_mass() const {return mass;}
    // Neutron Print
    void print_name()
    {
      cout<<"Neutron: [ID, mass, charge, spin, is anti] "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<"]"<<endl;
      cout<<"Expected energy Calorimeter: [EM_1: "<<energies_layers[0]<<", EM_2: "<<energies_layers[1]<<", HAD_1: "<<energies_layers[2]<<", HAD_2: "<<energies_layers[3]<<"]."<<endl;
      cout<<"Expected energy Tracker: [Inner: "<<energies_layers[4]<<", Outer: "<<energies_layers[5]<<", Strip: "<<energies_layers[6]<<"]"<<endl;
      cout<<"Expected energy Muon Chamber: [Inner: "<<energies_layers[7]<<", Outer: "<<energies_layers[8]<<"]"<<endl;
      cout<<"The true energy is: "<<true_energy<<endl;
    }
};
// END OF DIRECT PARTICLE CHILDREN
class Electron : public Lepton
{
  protected:
    vector<double> energies_layers{0,0,0,0,0,0,0,0,0};
  public:
    Electron()
    {
      mass = 0.511;
      charge = -1;
      spin = 0.5;
      true_energy = mass;

      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[2]; //TRACK1
      energies_layers[5] = divisions_true_energy[3]; //TRACK2
      energies_layers[6] = divisions_true_energy[4]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    Electron(bool is_anti_initial, double true_energy_initial) 
    {
      is_anti = is_anti_initial;
      if(is_anti == false)
      {
        charge = -1;
      }   else if(is_anti == true)
      {
        charge = 1;    
      }
      mass = 0.511;
      spin = 0.5;
      true_energy = mass;

      // Try and catch implementation...
      if(true_energy_initial >= mass)
      {
        true_energy = true_energy_initial;
        vector<double> divisions_true_energy{divide_random(true_energy, 5)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = 0; //HAD1
        energies_layers[3] = 0; //HAD2
        energies_layers[4] = divisions_true_energy[2]; //TRACK1
        energies_layers[5] = divisions_true_energy[3]; //TRACK2
        energies_layers[6] = divisions_true_energy[4]; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2
        return;
      }
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }

      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[2]; //TRACK1
      energies_layers[5] = divisions_true_energy[3]; //TRACK2
      energies_layers[6] = divisions_true_energy[4]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    ~Electron() {/*cout<<"Electron Destroyed"<<endl;*/}
    // Copy Constructor
    Electron(Electron const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
    }
    // Move Constructor
    Electron(Electron && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      } 
      // Now we clear
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
      energies_layers.clear(); // OJO AQUI
    }
    // Assignment operator
    Electron & operator=(Electron const & object)
    {
      if(&object == this) return *this;
      mass = object.mass;
      charge = object.charge;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
      return *this;
    }
    // Move Assignment Operator
    Electron & operator=(Electron && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      for(int i{0}; i < 9; i++)
      {
        std::swap(energies_layers[i], object.energies_layers[i]);
      }
      return *this;
    }

    // Setters
    void set_energies_in_layers(const double set_energy0, const double set_energy1, const double set_energy2, const double set_energy3)
    {
      if(set_energy0 + set_energy1 + set_energy2 + set_energy3 == true_energy and (set_energy0 > 0) and (set_energy1 > 0) and (set_energy2 > 0) and (set_energy3 > 0))
      {
        energies_layers[0] = set_energy0;
        energies_layers[1] = set_energy1;
        energies_layers[2] = set_energy2;
        energies_layers[3] = set_energy3;
      } else
      {
        cout<<"The values given does not sum up to the energy. The energy/c of the electron"<<endl;
        cout<<"will be divided equally between the four layers."<<endl;
        energies_layers[0] = true_energy/4;
        energies_layers[1] = true_energy/4;
        energies_layers[2] = true_energy/4;
        energies_layers[3] = true_energy/4;
      }
    }
    void set_true_energy(const double true_energy_set)
    {
      if(true_energy_set >= mass)
      {
      true_energy = true_energy_set;
      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[2]; //TRACK1
      energies_layers[5] = divisions_true_energy[3]; //TRACK2
      energies_layers[6] = divisions_true_energy[4]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
      return;
      }
      // Try and catch implementation...
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin>>new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }
      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[2]; //TRACK1
      energies_layers[5] = divisions_true_energy[3]; //TRACK2
      energies_layers[6] = divisions_true_energy[4]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
      
    }
    
    // Getters
    double energy_layer_n(const int n)
    {
      if(n == 0 or n == 1 or n == 2 or n == 3 or n == 4 or n == 5 or n == 6 or n == 7 or n == 8)
      {
        return energies_layers[n];
      } else
      {
        cout<<"Index out of bounds. Returning the content of the first layer."<<endl;
        return energies_layers[0];
      }
    }
    double get_mass() const {return mass;}
    int get_charge() const {return charge;}
    // Print Electron
    void print_name()
    {
      cout<<"Electron: [ID, mass, charge, spin, is anti] "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<"]"<<endl;
      cout<<"Expected energy Calorimeter: [EM_1: "<<energies_layers[0]<<", EM_2: "<<energies_layers[1]<<", HAD_1: "<<energies_layers[2]<<", HAD_2: "<<energies_layers[3]<<"]."<<endl;
      cout<<"Expected energy Tracker: [Inner: "<<energies_layers[4]<<", Outer: "<<energies_layers[5]<<", Strip: "<<energies_layers[6]<<"]"<<endl;
      cout<<"Expected energy Muon Chamber: [Inner: "<<energies_layers[7]<<", Outer: "<<energies_layers[8]<<"]"<<endl;
      cout<<"The true energy is: "<<true_energy<<endl;
    }

    bool get_decayed_yet() {cout<<"Electrons do not decay. Returning false."<<endl; return false;}
    void decay() {cout<<"Electrons do not decay"<<endl;}
};

class Muon : public Lepton
{
  protected:
    bool isolated{true};
    vector<double> energies_layers{0,0,0,0,0,0,0,0,0};

  public:
    Muon() 
    {
      charge = -1; mass = 105.66; true_energy = mass;

      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = 0; //EM1
      energies_layers[1] = 0; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[0]; //TRACK1
      energies_layers[5] = divisions_true_energy[1]; //TRACK2
      energies_layers[6] = divisions_true_energy[2]; //TRACK3
      energies_layers[7] = divisions_true_energy[3]; //MUON1
      energies_layers[8] = divisions_true_energy[4]; //MUON2
    }
    Muon(bool is_anti_initial, bool isolated_initial, double true_energy_initial)
    {
      mass = 105.66;
      spin = 0.5;
      isolated = isolated_initial;
      is_anti = is_anti_initial;
      if(is_anti == false)
      {
        charge = -1;
      } else if(is_anti == true)
      {
        charge = 1;    
      }
      true_energy = mass;

      // Try and catch implementation...
      if(true_energy_initial >= mass)
      {
        true_energy = true_energy_initial;
        vector<double> divisions_true_energy{divide_random(true_energy, 5)};
        energies_layers[0] = 0; //EM1
        energies_layers[1] = 0; //EM2
        energies_layers[2] = 0; //HAD1
        energies_layers[3] = 0; //HAD2
        energies_layers[4] = divisions_true_energy[0]; //TRACK1
        energies_layers[5] = divisions_true_energy[1]; //TRACK2
        energies_layers[6] = divisions_true_energy[2]; //TRACK3
        energies_layers[7] = divisions_true_energy[3]; //MUON1
        energies_layers[8] = divisions_true_energy[4]; //MUON2
        return;
      }
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }

      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = 0; //EM1
      energies_layers[1] = 0; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[0]; //TRACK1
      energies_layers[5] = divisions_true_energy[1]; //TRACK2
      energies_layers[6] = divisions_true_energy[2]; //TRACK3
      energies_layers[7] = divisions_true_energy[3]; //MUON1
      energies_layers[8] = divisions_true_energy[4]; //MUON2
    }
    ~Muon() {/*cout<<"Muon Destroyed"<<endl;*/}
    // Copy Constructor
    Muon(Muon const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      isolated = object.isolated;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
    }
    // Move Constructor
    Muon(Muon && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      isolated = object.isolated;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      } 
      // Now we clear
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
      isolated = false;
      energies_layers.clear(); // OJO AQUI
    }
    // Assignment operator
    Muon & operator=(Muon const & object)
    {
      if(&object == this) return *this;
      mass = object.mass;
      charge = object.charge;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      isolated = object.isolated;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
      return *this;
    }
    // Move Assignment Operator
    Muon & operator=(Muon && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      std::swap(isolated, object.isolated);
      for(int i{0}; i < 9; i++)
      {
        std::swap(energies_layers[i], object.energies_layers[i]);
      }
      return *this;
    }

    // Setters
    void set_isolation(const bool isolation_set) {isolated = isolation_set;}
    void set_true_energy(const double true_energy_set)
    {
      if(true_energy_set >= mass)
      {
      true_energy = true_energy_set;
      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = 0; //EM1
      energies_layers[1] = 0; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[0]; //TRACK1
      energies_layers[5] = divisions_true_energy[1]; //TRACK2
      energies_layers[6] = divisions_true_energy[2]; //TRACK3
      energies_layers[7] = divisions_true_energy[3]; //MUON1
      energies_layers[8] = divisions_true_energy[4]; //MUON2
      return;
      }
      // Try and catch implementation...
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }
      vector<double> divisions_true_energy{divide_random(true_energy, 5)};
      energies_layers[0] = 0; //EM1
      energies_layers[1] = 0; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = divisions_true_energy[0]; //TRACK1
      energies_layers[5] = divisions_true_energy[1]; //TRACK2
      energies_layers[6] = divisions_true_energy[2]; //TRACK3
      energies_layers[7] = divisions_true_energy[3]; //MUON1
      energies_layers[8] = divisions_true_energy[4]; //MUON2
    }

    // Getters
    bool get_isolation() const {return isolated;}
    int get_charge() const {return charge;}
    double get_mass() const {return mass;}
    vector<int> get_muon_seeds() const {return muon_seed;}
    double energy_layer_n(const int n)
    {
      if(n == 0 or n == 1 or n == 2 or n == 3 or n == 4 or n == 5 or n== 6 or n == 7 or n == 8)
      {
        return energies_layers[n];
      } else
      {
        cout<<"Index out of bounds. Returning the content of the first layer."<<endl;
        return energies_layers[0];
      }
    }

    // Print Muon
    void print_name()
    {
      cout<<"Muon: [ID, mass, charge, spin, is anti] "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<"]"<<endl;
      cout<<"Expected energy Calorimeter: [EM_1: "<<energies_layers[0]<<", EM_2: "<<energies_layers[1]<<", HAD_1: "<<energies_layers[2]<<", HAD_2: "<<energies_layers[3]<<"]."<<endl;
      cout<<"Expected energy Tracker: [Inner: "<<energies_layers[4]<<", Outer: "<<energies_layers[5]<<", Strip: "<<energies_layers[6]<<"]"<<endl;
      cout<<"Expected energy Muon Chamber: [Inner: "<<energies_layers[7]<<", Outer: "<<energies_layers[8]<<"]"<<endl;
      cout<<"Isolated Muon: "<<isolated<<endl;
      cout<<"The true energy is: "<<true_energy<<endl;
    }

    bool get_decayed_yet() {cout<<"Muons do not decay. Returning false."<<endl; return false;}
    void decay() {cout<<"Muons do not decay"<<endl;}
};

class Neutrino : public Lepton
{
  protected:
    char flavour{'e'}; // Can take value e, m, t 
    bool has_interacted{false};

  public:
    Neutrino() {mass = 0; charge = 0; spin = 0.5;}
    Neutrino(char flavour_initial, bool has_interacted_initial, bool is_anti_initial)
    {
      if(flavour_initial != 'e' and flavour_initial != 'm' and flavour_initial != 't')
      {
        cout<<"The Neutrino type must be electron (e), muon (m) or tau (t)."<<endl;
        cin>>flavour_initial;
        while(cin.fail() or (flavour_initial != 'e' and flavour_initial != 'm' and flavour_initial != 't'))
        {
          cout<<"Invalid input, try again: ";
          cin.clear();
          cin.ignore();
          cin>>flavour_initial;
        }
        flavour = flavour_initial;
      } else
      {
        flavour = flavour_initial;
      }
      mass = 0;
      charge = 0;
      is_anti = is_anti_initial;
      spin = 0.5;
    }
    ~Neutrino() {cout<<"Neutrino Destroyed"<<endl;}
    // Copy Constructor
    Neutrino(Neutrino const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      flavour = object.flavour;
      has_interacted = object.has_interacted;
    }
    // Move Constructor
    Neutrino(Neutrino && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      flavour = object.flavour;
      has_interacted = object.has_interacted;
      // Now we clear
      charge = 0;
      flavour = 'e';
      has_interacted = false;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
    }
    // Assignmet Operator
    Neutrino & operator=(Neutrino const & object)
    {
      if(&object == this) return *this;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      flavour = object.flavour;
      has_interacted = object.has_interacted;
      return *this;
    }
    // Move Assignment Operator
    Neutrino & operator=(Neutrino && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      std::swap(has_interacted, object.has_interacted);
      std::swap(flavour, object.flavour);
      return *this;
    }

    // Setters
    void set_flavour(const char flavour_set)
    {
      if (flavour_set != 'e' and flavour_set != 'm' and flavour_set != 't')
      {
        cout<<"Flavour selected does not exist. The deafult is electron (e)."<<endl;
        flavour = 'e';
      } else
      {
        flavour =  flavour_set;
      }
    }
    void set_interaction(const bool interaction_set) {has_interacted =  interaction_set;}

    // Getters
    bool get_interacted() const {return has_interacted;}
    char get_flavour() const {return flavour;}
    int get_charge() const {return charge;}
    double get_mass() const {return mass;}
    // Print Neutrino
    void print_name()
      {
        cout<<"Neutrino: [ID, mass, charge, spin, antiparticle, flavour, Interacted] = "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<", "<<flavour<<", "<<has_interacted<<"]. ";
        cout<<"The True energy is: "<<true_energy<<endl;
      }

    bool get_decayed_yet() {cout<<"Neutrinos do not decay. Returning false."<<endl; return false;}
    void decay() {cout<<"Neutrinos do not decay"<<endl;}
};

class Tau : public Lepton
{
  protected:
    string decay_type{"Undefined"}; // set to hadronic or leptonic
    vector<std::shared_ptr<Particle>> leptonic_decay_products; // leptonic decay products
    vector<double> energies_layers{0,0,0,0,0,0,0,0,0};
    int decay_random{rand() % 10 + 1};
    bool decayed_yet{false};

  public:
    Tau() 
    {
      mass = 1776.86; 
      spin = 0.5;
      true_energy = mass;
      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2

      // Making leptonic or hadronic randomly
      if(rand() % 10 + 1 > likelihood_leptonic_decay)
      {
        decay_type = "leptonic";
      } else
      {
        decay_type = "hadronic";
      }
      if(is_anti == false)
      {
        charge = -1; 
      } else
      {
        charge = 1;
      }
    };
    Tau(bool is_anti_initial, string decay_type_initial, double true_energy_initial)
    {
      if(decay_type_initial != "leptonic" and decay_type_initial != "hadronic")
      {
        cout<<"The decay type of Tau should be Leptonic or Hadronic."<<endl;
        cin>> decay_type_initial;
        while(cin.fail() or (decay_type_initial != "leptonic" and decay_type_initial != "hadronic"))
        {
          cout<<"Invalid input, try again: ";
          cin.clear();
          cin.ignore();
          cin>>decay_type_initial;
        }
        decay_type = decay_type_initial;
      } else
      {
        decay_type = decay_type_initial;
      }
      is_anti = is_anti_initial;
      mass = 1776.86;
      spin = 0.5;
      true_energy = mass;
      if(is_anti == false)
      {
        charge = -1;
      } else 
      {
        charge = 1;
      }

      // Try and catch implementation...
      if(true_energy_initial >= mass)
      {
        true_energy = true_energy_initial;
        vector<double> divisions_true_energy{divide_random(true_energy, 7)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = divisions_true_energy[2]; //HAD1
        energies_layers[3] = divisions_true_energy[3]; //HAD2
        energies_layers[4] = divisions_true_energy[4]; //TRACK1
        energies_layers[5] = divisions_true_energy[5]; //TRACK2
        energies_layers[6] = divisions_true_energy[6]; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2
        return;
      }
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }

      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    ~Tau() {cout<<"Tau Destroyed"<<endl;}
    // Copy Constructor
    Tau(Tau const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      decay_type = object.decay_type;
      decay_random = object.decay_random;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
      if(decayed_yet == false and object.decayed_yet == true)
      {
        leptonic_decay_products.push_back(object.leptonic_decay_products[0]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[1]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[2]);
        decayed_yet = true;
      } else if(decayed_yet == true and object.decayed_yet == true)
      {
        leptonic_decay_products[0] = object.leptonic_decay_products[0];
        leptonic_decay_products[1] = object.leptonic_decay_products[1];
        leptonic_decay_products[2] = object.leptonic_decay_products[2];
      } else if(decayed_yet == false and object.decayed_yet == false)
      {
        decayed_yet = false;
      } else if(decayed_yet == true and object.decayed_yet == false)
      {
        decayed_yet = false;
        leptonic_decay_products.clear();
      }
    }
    // Move Constructor
    Tau(Tau && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      decay_type = object.decay_type;
      decay_random = object.decay_random;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      } 
      if(decayed_yet == false and object.decayed_yet == true)
      {
        leptonic_decay_products.push_back(object.leptonic_decay_products[0]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[1]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[2]);
        decayed_yet = true;
      } else if(decayed_yet == true and object.decayed_yet == true)
      {
        leptonic_decay_products[0] = object.leptonic_decay_products[0];
        leptonic_decay_products[1] = object.leptonic_decay_products[1];
        leptonic_decay_products[2] = object.leptonic_decay_products[2];
      } else if (decayed_yet == false and object.decayed_yet == false)
      {
        decayed_yet = false;
      } else if(decayed_yet == true and object.decayed_yet == false)
      {
        decayed_yet = false;
      }

      // Now we clear
      decay_random = 0;
      decay_type = "undefined";
      object.leptonic_decay_products.clear();
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
      energies_layers.clear(); // OJO AQUI
    }
    // Assignment Operator
    Tau & operator=(Tau const & object)
    {
      if(&object == this) return *this;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      decay_type = object.decay_type;
      decay_random = object.decay_random;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
      if(decayed_yet == false and object.decayed_yet == true)
      {
        leptonic_decay_products.push_back(object.leptonic_decay_products[0]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[1]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[2]);
        decayed_yet = true;
      } else if(decayed_yet == true and object.decayed_yet == true)
      {
        leptonic_decay_products[0] = object.leptonic_decay_products[0];
        leptonic_decay_products[1] = object.leptonic_decay_products[1];
        leptonic_decay_products[2] = object.leptonic_decay_products[2];
      } else if (decayed_yet == false and object.decayed_yet == false)
      {
        decayed_yet = false;
      }
      return *this;
    }
    // Move Assignment Operator
    Tau & operator=(Tau && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      std::swap(decay_type, object.decay_type);
      std::swap(decay_random, object.decay_random);
      for(int i{0}; i < 9; i++)
      {
        std::swap(energies_layers[i], object.energies_layers[i]);
      }
      if(decayed_yet == false and object.decayed_yet == true)
      {
        leptonic_decay_products.push_back(object.leptonic_decay_products[0]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[1]);
        leptonic_decay_products.push_back(object.leptonic_decay_products[2]);
        decayed_yet = true;
      } else if(decayed_yet == true and object.decayed_yet == true)
      {
        leptonic_decay_products[0] = object.leptonic_decay_products[0];
        leptonic_decay_products[1] = object.leptonic_decay_products[1];
        leptonic_decay_products[2] = object.leptonic_decay_products[2];
      } else if (decayed_yet == false and object.decayed_yet == false)
      {
        decayed_yet = false;
      } else if(decayed_yet == true and object.decayed_yet == false)
      {
        decayed_yet = false;
        leptonic_decay_products.clear();
      }
      return *this;
    }

    // Setters
    void set_decay(const string decay_set)
    {
      if(decay_set == "leptonic" or decay_set == "hadronic")
      {
        decay_type = decay_set;
      } else
      {
        cout<<"The decay type should be either leptonic or hadronic. The default value is leptonic."<<endl;
        decay_type = "leptonic";
      }
    }
    void set_true_energy(const double true_energy_set)
    {
      if(true_energy_set >= mass)
      {
        true_energy = true_energy_set;
        vector<double> divisions_true_energy{divide_random(true_energy, 7)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = divisions_true_energy[2]; //HAD1
        energies_layers[3] = divisions_true_energy[3]; //HAD2
        energies_layers[4] = divisions_true_energy[4]; //TRACK1
        energies_layers[5] = divisions_true_energy[5]; //TRACK2
        energies_layers[6] = divisions_true_energy[6]; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2  
        return;
      }
      // Try and catch implementation...
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }
      vector<double> divisions_true_energy{divide_random(true_energy, 7)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = divisions_true_energy[2]; //HAD1
      energies_layers[3] = divisions_true_energy[3]; //HAD2
      energies_layers[4] = divisions_true_energy[4]; //TRACK1
      energies_layers[5] = divisions_true_energy[5]; //TRACK2
      energies_layers[6] = divisions_true_energy[6]; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    
    // Getters
    string get_decay_type() const {return decay_type;}
    vector<std::shared_ptr<Particle>> get_decay_products() {return leptonic_decay_products;}
    int get_charge() const {return charge;}
    double get_mass() const {return mass;}
    bool get_decayed_yet() const override {return decayed_yet;}
    double energy_layer_n(const int n)
    {
      if(n == 0 or n == 1 or n == 2 or n == 3 or n == 4 or n == 5 or n == 6 or n == 7 or n == 8)
      {
        return energies_layers[n];
      } else
      {
        cout<<"Index out of bounds. Returning the content of the first layer."<<endl;
        return energies_layers[0];
      }
    }

   // Decay Function
   // Once the Tau decays then its leptonic_decay_products vector is filled.
   void decay()
   {
      if(decayed_yet == true)
      {
        return;
      }
      decayed_yet = true;
      if(decay_type == "hadronic")
      {
        cout<<"The tau decays hadronically..."<<endl;
        return;
      }
      // Preparatin to distribute the energy of the Tau randomly between decay products.
      vector<double> divs_energy{divide_random(true_energy, 3)};

      if(is_anti == false)
      {
        if(decay_random < 5)
        {
          leptonic_decay_products.push_back(std::make_shared<Electron>(0, divs_energy[0]));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('e', 0, 1));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('t', 0, 0));
          leptonic_decay_products[0]->set_true_energy(divs_energy[0]);
          leptonic_decay_products[1]->set_true_energy(divs_energy[1]);
          leptonic_decay_products[2]->set_true_energy(divs_energy[2]);
          cout<<"Tau Decayed."<<endl;
        }   else
        {
          leptonic_decay_products.push_back(std::make_shared<Muon>(0, 1, divs_energy[0]));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('m', 0, 1));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('t', 0, 0));
          leptonic_decay_products[0]->set_true_energy(divs_energy[0]);
          leptonic_decay_products[1]->set_true_energy(divs_energy[1]);
          leptonic_decay_products[2]->set_true_energy(divs_energy[2]);
          cout<<"Tau Decayed."<<endl;
        }
      } else
      {
        if(decay_random < 5)
        {
          leptonic_decay_products.push_back(std::make_shared<Electron>(1, divs_energy[0]));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('e', 0, 0));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('t', 0, 1));
          leptonic_decay_products[0]->set_true_energy(divs_energy[0]);
          leptonic_decay_products[1]->set_true_energy(divs_energy[1]);
          leptonic_decay_products[2]->set_true_energy(divs_energy[2]);
          cout<<"Tau Decayed."<<endl;
        }   else
        {
          leptonic_decay_products.push_back(std::make_shared<Muon>(1, 1, divs_energy[0]));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('m', 0, 0));
          leptonic_decay_products.push_back(std::make_shared<Neutrino>('t', 0, 1));
          leptonic_decay_products[0]->set_true_energy(divs_energy[0]);
          leptonic_decay_products[1]->set_true_energy(divs_energy[1]);
          leptonic_decay_products[2]->set_true_energy(divs_energy[2]);
          cout<<"Tau Decayed."<<endl;
        }
      }
   }

   // Print Tau
   void print_name()
   {
      cout<<"Tau: [ID, mass, charge, spin, is anti] "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<", "<<"Anti: "<<is_anti<<"]"<<endl;
      cout<<"Decay Type: "<<decay_type<<endl;
      cout<<"Expected energy Calorimeter: [EM_1: "<<energies_layers[0]<<", EM_2: "<<energies_layers[1]<<", HAD_1: "<<energies_layers[2]<<", HAD_2: "<<energies_layers[3]<<"]."<<endl;
      cout<<"Expected energy Tracker: [Inner: "<<energies_layers[4]<<", Outer: "<<energies_layers[5]<<", Strip: "<<energies_layers[6]<<"]"<<endl;
      cout<<"Expected energy Muon Chamber: [Inner: "<<energies_layers[7]<<", Outer: "<<energies_layers[8]<<"]"<<endl;
      cout<<"The true energy is: "<<true_energy<<endl;
   }
};

// END OF DIRECT LEPTON CHILDREN
class Photon : public Boson
{
  protected:
    vector<double> energies_layers{0,0,0,0,0,0,0,0,0};
  public:
    Photon()
    {
      mass = 0;
      charge = 0;
      spin = 1;
      true_energy = (rand() % 10 + 1) * 100; // From 100 to 1000.
      vector<double> divisions_true_energy{divide_random(true_energy, 2)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = 0; //TRACK1
      energies_layers[5] = 0; //TRACK2
      energies_layers[6] = 0; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    Photon(double true_energy_initial)
    {
      mass = 0;
      charge = 0;
      spin = 1;
      
      // Try and catch implementation...
      if(true_energy_initial >= mass)
      {
        true_energy = true_energy_initial;
        vector<double> divisions_true_energy{divide_random(true_energy, 2)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = 0; //HAD1
        energies_layers[3] = 0; //HAD2
        energies_layers[4] = 0; //TRACK1
        energies_layers[5] = 0; //TRACK2
        energies_layers[6] = 0; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2
        return;
      }
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }

      vector<double> divisions_true_energy{divide_random(true_energy, 2)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = 0; //TRACK1
      energies_layers[5] = 0; //TRACK2
      energies_layers[6] = 0; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    ~Photon() {cout<<"Photon Destroyed"<<endl;}
    // Copy Constructor
    Photon(Photon const & object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
    }
    // Move Constructor
    Photon(Photon && object)
    {
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      charge = object.charge;
      mass = object.mass;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      } 
      // Now we clear
      charge = 0;
      mass = 0; 
      identifier = -1; // Unconsidered
      true_energy = 0;
      spin = 0;
      is_anti = false;
      tracker_seeds.clear();
      muon_seed.clear();
      calorimeter_seed = 0;
      energies_layers.clear(); // OJO AQUI
    }
    // Assignment operator
    Photon & operator=(Photon const & object)
    {
      if(&object == this) return *this;
      mass = object.mass;
      charge = object.charge;
      identifier = object.identifier;
      true_energy = object.true_energy;
      spin = object.spin;
      is_anti = object.is_anti;
      for(int i{0}; i < 9; i++)
      {
        energies_layers[i] = object.energies_layers[i];
      }
      return *this;
    }
    // Move Assignment Operator
    Photon & operator=(Photon && object)
    {
      std::swap(identifier, object.identifier);
      std::swap(true_energy, object.true_energy);
      std::swap(spin, object.spin);
      std::swap(is_anti, object.is_anti);
      std::swap(mass, object.mass);
      std::swap(charge, object.charge);
      for(int i{0}; i < 9; i++)
      {
        std::swap(energies_layers[i], object.energies_layers[i]);
      }
      return *this;
    }

    // Setters
    void set_energies_in_layers(const double set_energy0, const double set_energy1, const double set_energy2, const double set_energy3)
    {
      if(set_energy0 + set_energy1 + set_energy2 + set_energy3 == true_energy and (set_energy0 > 0) and (set_energy1 > 0) and (set_energy2 > 0) and (set_energy3 > 0))
      {
        energies_layers[0] = set_energy0;
        energies_layers[1] = set_energy1;
        energies_layers[2] = set_energy2;
        energies_layers[3] = set_energy3;
      } else
      {
        cout<<"The values given does not sum up to the energy or one of them is negative. The energy/c of the electron"<<endl;
        cout<<"will be divided equally between the four layers."<<endl;
        energies_layers[0] = true_energy/4;
        energies_layers[1] = true_energy/4;
        energies_layers[2] = true_energy/4;
        energies_layers[3] = true_energy/4;
      }
    }
    void set_mass(const double mass_set)
    {
      cout<<"Photons are massless"<<endl;
      mass = 0;
    }
    void set_charge()
    {
      cout<<"Photons do not carry charge"<<endl;
      charge = 0;
    }
    void set_true_energy(const double true_energy_set)
    {
      if(true_energy_set >= mass)
      {
        true_energy = true_energy_set;
        vector<double> divisions_true_energy{divide_random(true_energy, 2)};
        energies_layers[0] = divisions_true_energy[0]; //EM1
        energies_layers[1] = divisions_true_energy[1]; //EM2
        energies_layers[2] = 0; //HAD1
        energies_layers[3] = 0; //HAD2
        energies_layers[4] = 0; //TRACK1
        energies_layers[5] = 0; //TRACK2
        energies_layers[6] = 0; //TRACK3
        energies_layers[7] = 0; //MUON1
        energies_layers[8] = 0; //MUON2
        return;
      }
      // Try and catch implementation...
      double new_energy;
      bool valid_input{false};
      while (valid_input != true)
      {
        try
        {
          cout<<"Try introducing another energy:";
          cin >> new_energy;
          if(new_energy <= mass)
          {
            throw std::invalid_argument(" ----> Energy must be greater than the invariant mass of the particle <---- ");
          } else
          {
            valid_input = true;
            true_energy = new_energy;
          }
        }
        catch(const std::exception& e)
        {
          std::cerr<<"Error: "<<e.what()<<endl;
          cin.clear();
          cin.ignore();
        } 
      }
      vector<double> divisions_true_energy{divide_random(true_energy, 2)};
      energies_layers[0] = divisions_true_energy[0]; //EM1
      energies_layers[1] = divisions_true_energy[1]; //EM2
      energies_layers[2] = 0; //HAD1
      energies_layers[3] = 0; //HAD2
      energies_layers[4] = 0; //TRACK1
      energies_layers[5] = 0; //TRACK2
      energies_layers[6] = 0; //TRACK3
      energies_layers[7] = 0; //MUON1
      energies_layers[8] = 0; //MUON2
    }
    // Getters
    double energy_layer_n(const int n)
    {
      if(n == 0 or n == 1 or n == 2 or n == 3 or n == 4 or n == 5 or n == 6 or n == 7 or n == 8)
      {
        return energies_layers[n];
      } else
      {
        cout<<"Index out of bounds. Returning the content of the first layer."<<endl;
        return energies_layers[0];
      }
    }
    double get_mass() const { return mass;}
    int get_charge() {return charge;}
    // Photon Print
    void print_name()
    {
    cout<<"Photon: [ID, mass, charge, spin] "<<"["<<identifier<<", "<<mass<<" Mev/c^2"<<", "<<charge<<"e"<<", "<<spin<<"]"<<endl;
    cout<<"Expected energy Calorimeter: [EM_1: "<<energies_layers[0]<<", EM_2: "<<energies_layers[1]<<", HAD_1: "<<energies_layers[2]<<", HAD_2: "<<energies_layers[3]<<"]."<<endl;
    cout<<"Expected energy Tracker: [Inner: "<<energies_layers[4]<<", Outer: "<<energies_layers[5]<<", Strip: "<<energies_layers[6]<<"]"<<endl;
    cout<<"Expected energy Muon Chamber: [Inner: "<<energies_layers[7]<<", Outer: "<<energies_layers[8]<<"]"<<endl;
    cout<<"The true energy is: "<<true_energy<<endl;
    }
};


// DETECTOR CLASSES
class Tracker
{
  protected:
    bool is_on{true};
    vector<std::shared_ptr<Particle>> INNER_PIXEL_LAYER_particles; // We store the particles detected by each layer in this vectors 
    vector<std::shared_ptr<Particle>> OUTER_PIXEL_LAYER_particles; // If one particle appears in two of these then the particle is considered to be detected.
    vector<std::shared_ptr<Particle>> STRIP_LAYER_particles;
    vector<double> INNER_PIXEL_LAYER_energies; // ESTAS LINEAS PUEDEN SER INNECESARIAS SI SOLO EL
    vector<double> OUTER_PIXEL_LAYER_energies; // CALORIMETRO DETECTA LAS ENERGIAS.
    vector<double> STRIP_LAYER_energies;
    vector<std::shared_ptr<Particle>> particles_detected;

  public:
    Tracker() = default;
    Tracker(bool is_on_initial) {is_on = is_on_initial;}
    ~Tracker() {cout<<"Tracker Destroyed"<<endl;}
    // Copy Constructor
    Tracker(Tracker const & object)
    {
      is_on = object.is_on;
      // Make sure all lists are empty
      INNER_PIXEL_LAYER_particles.clear();
      OUTER_PIXEL_LAYER_particles.clear();
      STRIP_LAYER_particles.clear();
      INNER_PIXEL_LAYER_energies.clear();
      OUTER_PIXEL_LAYER_energies.clear();
      STRIP_LAYER_energies.clear();
      particles_detected.clear();

      // Fill with content of object
      for(std::shared_ptr<Particle> particle : object.INNER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.STRIP_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_PIXEL_LAYER_energies)
      {
        INNER_PIXEL_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_PIXEL_LAYER_energies)
      {
        OUTER_PIXEL_LAYER_energies.push_back(energy );
      }
      for(double energy : object.STRIP_LAYER_energies)
      {
        STRIP_LAYER_energies .push_back(energy );
      }
    }
    // Move Constructor
    Tracker(Tracker && object)
    {
      is_on = object.is_on;
      // Make sure all lists are empty
      INNER_PIXEL_LAYER_particles.clear();
      OUTER_PIXEL_LAYER_particles.clear();
      STRIP_LAYER_particles.clear();
      INNER_PIXEL_LAYER_energies.clear();
      OUTER_PIXEL_LAYER_energies.clear();
      STRIP_LAYER_energies.clear();
      particles_detected.clear();

      // Fill with content of object
      for(std::shared_ptr<Particle> particle : object.INNER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.STRIP_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_PIXEL_LAYER_energies)
      {
        INNER_PIXEL_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_PIXEL_LAYER_energies)
      {
        OUTER_PIXEL_LAYER_energies.push_back(energy );
      }
      for(double energy : object.STRIP_LAYER_energies)
      {
        STRIP_LAYER_energies .push_back(energy );
      }
      // Now we clear the other object
      object.is_on = false;
      object.INNER_PIXEL_LAYER_particles.clear();
      object.OUTER_PIXEL_LAYER_particles.clear();
      object.STRIP_LAYER_particles.clear();
      object.INNER_PIXEL_LAYER_energies.clear();
      object.OUTER_PIXEL_LAYER_energies.clear();
      object.STRIP_LAYER_energies.clear();
      object.particles_detected.clear();
    }
    // Assignment Operator
    Tracker & operator=(Tracker const & object)
    {
      if(&object == this) return *this;
      is_on = object.is_on;
      // Make sure all lists are empty
      INNER_PIXEL_LAYER_particles.clear();
      OUTER_PIXEL_LAYER_particles.clear();
      STRIP_LAYER_particles.clear();
      INNER_PIXEL_LAYER_energies.clear();
      OUTER_PIXEL_LAYER_energies.clear();
      STRIP_LAYER_energies.clear();
      particles_detected.clear();

      // Fill with content of object
      for(std::shared_ptr<Particle> particle : object.INNER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.STRIP_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_PIXEL_LAYER_energies)
      {
        INNER_PIXEL_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_PIXEL_LAYER_energies)
      {
        OUTER_PIXEL_LAYER_energies.push_back(energy );
      }
      for(double energy : object.STRIP_LAYER_energies)
      {
        STRIP_LAYER_energies .push_back(energy );
      }
      return *this;
    }
    // Move Assignment Operator
    Tracker & operator=(Tracker && object)
    {
      std::swap(is_on, object.is_on);
      // Make sure all lists are empty
      INNER_PIXEL_LAYER_particles.clear();
      OUTER_PIXEL_LAYER_particles.clear();
      STRIP_LAYER_particles.clear();
      INNER_PIXEL_LAYER_energies.clear();
      OUTER_PIXEL_LAYER_energies.clear();
      STRIP_LAYER_energies.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.INNER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_PIXEL_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.STRIP_LAYER_particles)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_PIXEL_LAYER_energies)
      {
        INNER_PIXEL_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_PIXEL_LAYER_energies)
      {
        OUTER_PIXEL_LAYER_energies.push_back(energy );
      }
      for(double energy : object.STRIP_LAYER_energies)
      {
        STRIP_LAYER_energies .push_back(energy );
      }
      // Now we clear the other object
      object.is_on = false;
      object.INNER_PIXEL_LAYER_particles.clear();
      object.OUTER_PIXEL_LAYER_particles.clear();
      object.STRIP_LAYER_particles.clear();
      object.INNER_PIXEL_LAYER_energies.clear();
      object.OUTER_PIXEL_LAYER_energies.clear();
      object.STRIP_LAYER_energies.clear();
      object.particles_detected.clear();
      return *this;
    }

    // Setters
    void set_is_on(bool is_on_set) {is_on = is_on_set;}

    // Getters
    vector<std::shared_ptr<Particle>> & get_particles_INNER_PIXEL_LAYER() {return INNER_PIXEL_LAYER_particles;} // ESTO ESTA PENDIENTE AUN
    vector<std::shared_ptr<Particle>> & get_particles_OUTER_PIXEL_LAYER() {return OUTER_PIXEL_LAYER_particles;} // ESTO ESTA PENDIENTE AUN
    vector<std::shared_ptr<Particle>> & get_particles_STRIP_LAYER() {return STRIP_LAYER_particles;} // ESTO ESTA PENDIENTE AUN
    vector<double> get_INNER_PIXEL_LAYER_energies() {return INNER_PIXEL_LAYER_energies;}
    vector<double> get_OUTER_PIXEL_LAYER_energies() {return OUTER_PIXEL_LAYER_energies;}
    vector<double> get_STRIP_LAYER_energies() {return STRIP_LAYER_energies;}
    vector<std::shared_ptr<Particle>> & get_particles_detected() {return particles_detected;}

    // Make default
    void clear_tracker()
    {
      INNER_PIXEL_LAYER_energies.clear();
      OUTER_PIXEL_LAYER_energies.clear();
      STRIP_LAYER_energies.clear();
      INNER_PIXEL_LAYER_particles.clear();
      OUTER_PIXEL_LAYER_particles.clear();
      STRIP_LAYER_particles.clear();
      particles_detected.clear();
    }
    // Print Tracker
    void print_name()
    {
      cout<<"Tracker: [total particles detected, INNER_PIXEL_LAYER, OUTER_PIXEL_LAYER, STRIP_LAYER] = ["<<particles_detected.size()<<", "<<INNER_PIXEL_LAYER_particles.size()<<", "<<OUTER_PIXEL_LAYER_particles.size()<<", "<<STRIP_LAYER_particles.size()<<"]"<<endl;
    }

    // Detection Process
    // -----------------
    bool detect(std::shared_ptr<Particle> & particle)
    {
      if((particle->get_charge() == 0 and particle->get_tracker_seed()[0] > 10 - tracker_efficiency) or particle->get_mass() == 0) // Neutral particles do not interact with the tracker.
      {
        // Change tracker random numbers: 
        particle->change_seeds_tracker(); // Useful for multidetection process
        return false;
      }

      int detected_sum{0}; // To see if the particle was detected by more than just one layer
      // LAYER DETECTION
      if(particle->get_tracker_seed()[0] > 10 - tracker_efficiency)
      {
        INNER_PIXEL_LAYER_particles.push_back(particle);
        INNER_PIXEL_LAYER_energies.push_back(particle->energy_layer_n(4));
        detected_sum = detected_sum + 1;
      }
      if(particle->get_tracker_seed()[1] > 10 - tracker_efficiency)
      {
        OUTER_PIXEL_LAYER_particles.push_back(particle);
        INNER_PIXEL_LAYER_energies.push_back(particle->energy_layer_n(5));
        detected_sum = detected_sum + 1;
      } 
      if(particle->get_tracker_seed()[2] > 10 - tracker_efficiency)
      {
        STRIP_LAYER_particles.push_back(particle);
        INNER_PIXEL_LAYER_energies.push_back(particle->energy_layer_n(6));
        detected_sum = detected_sum + 1;
      }
      particle->change_seeds_tracker();

      // At this point the particle is detected or not by the layers. It will return true for detection and false otherwise
      if(detected_sum >= 2)
      {
        particles_detected.push_back(particle);
        return true;
      } else
      {
        return false;
      }
    }
};

class Calorimeter
{
  protected:
    bool is_on{true};
    vector<unsigned int> EM_1_particles_ID; 
    vector<unsigned int> EM_2_particles_ID; 
    vector<unsigned int> HAD_1_particles_ID; 
    vector<unsigned int> HAD_2_particles_ID;
    vector<double> EM_1_energies;
    vector<double> EM_2_energies;
    vector<double> HAD_1_energies;
    vector<double> HAD_2_energies;
    vector<std::shared_ptr<Particle>> particles_detected;

  public:
    Calorimeter() = default;
    Calorimeter(bool is_on_initial) { is_on = is_on_initial;}
    ~Calorimeter() {cout<<"Calorimeter Destroyed"<<endl;}
    // Copy Constructor 
    Calorimeter(Calorimeter const & object)
    {
      is_on = object.is_on;
      // Make sure all vectors are empty
      EM_1_particles_ID.clear();
      EM_2_particles_ID.clear();
      HAD_1_particles_ID.clear();
      HAD_2_particles_ID.clear();
      EM_1_energies.clear();
      EM_2_energies.clear();
      HAD_1_energies.clear();
      HAD_2_energies.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(unsigned int id : object.EM_1_particles_ID)
      {
        EM_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.EM_2_particles_ID)
      {
        EM_2_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_1_particles_ID)
      {
        HAD_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_2_particles_ID)
      {
        HAD_2_particles_ID.push_back(id);
      }
      for(double energy : object.EM_1_energies)
      {
        EM_1_energies.push_back(energy);
      }
      for(double energy : object.EM_2_energies)
      {
        EM_2_energies.push_back(energy);
      }
      for(double energy : object.HAD_1_energies)
      {
        HAD_1_energies.push_back(energy);
      }
      for(double energy : object.HAD_2_energies)
      {
        HAD_2_energies.push_back(energy);
      }
    }

    // Move Constructor
    Calorimeter(Calorimeter && object)
    {
      is_on = object.is_on;
      // Make sure all vectors are empty
      EM_1_particles_ID.clear();
      EM_2_particles_ID.clear();
      HAD_1_particles_ID.clear();
      HAD_2_particles_ID.clear();
      EM_1_energies.clear();
      EM_2_energies.clear();
      HAD_1_energies.clear();
      HAD_2_energies.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(unsigned int id : object.EM_1_particles_ID)
      {
        EM_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.EM_2_particles_ID)
      {
        EM_2_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_1_particles_ID)
      {
        HAD_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_2_particles_ID)
      {
        HAD_2_particles_ID.push_back(id);
      }
      for(double energy : object.EM_1_energies)
      {
        EM_1_energies.push_back(energy);
      }
      for(double energy : object.EM_2_energies)
      {
        EM_2_energies.push_back(energy);
      }
      for(double energy : object.HAD_1_energies)
      {
        HAD_1_energies.push_back(energy);
      }
      for(double energy : object.HAD_2_energies)
      {
        HAD_2_energies.push_back(energy);
      }
      object.is_on = false;
      object.EM_1_particles_ID.clear();
      object.EM_2_particles_ID.clear();
      object.HAD_1_particles_ID.clear();
      object.HAD_2_particles_ID.clear();
      object.EM_1_energies.clear();
      object.EM_2_energies.clear();
      object.HAD_1_energies.clear();
      object.HAD_2_energies.clear();
      object.particles_detected.clear();
    }

    // Assignment Operator
    Calorimeter & operator=(Calorimeter const & object)
    {
      if(&object == this) return *this;
      is_on = object.is_on;
      // Make sure all vectors are empty
      EM_1_particles_ID.clear();
      EM_2_particles_ID.clear();
      HAD_1_particles_ID.clear();
      HAD_2_particles_ID.clear();
      EM_1_energies.clear();
      EM_2_energies.clear();
      HAD_1_energies.clear();
      HAD_2_energies.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(unsigned int id : object.EM_1_particles_ID)
      {
        EM_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.EM_2_particles_ID)
      {
        EM_2_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_1_particles_ID)
      {
        HAD_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_2_particles_ID)
      {
        HAD_2_particles_ID.push_back(id);
      }
      for(double energy : object.EM_1_energies)
      {
        EM_1_energies.push_back(energy);
      }
      for(double energy : object.EM_2_energies)
      {
        EM_2_energies.push_back(energy);
      }
      for(double energy : object.HAD_1_energies)
      {
        HAD_1_energies.push_back(energy);
      }
      for(double energy : object.HAD_2_energies)
      {
        HAD_2_energies.push_back(energy);
      }
      return *this;
    }

    // Move Assignment Operator
    Calorimeter & operator=(Calorimeter && object)
    {
      std::swap(is_on, object.is_on);
      EM_1_particles_ID.clear();
      EM_2_particles_ID.clear();
      HAD_1_particles_ID.clear();
      HAD_2_particles_ID.clear();
      EM_1_energies.clear();
      EM_2_energies.clear();
      HAD_1_energies.clear();
      HAD_2_energies.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(unsigned int id : object.EM_1_particles_ID)
      {
        EM_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.EM_2_particles_ID)
      {
        EM_2_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_1_particles_ID)
      {
        HAD_1_particles_ID.push_back(id);
      }
      for(unsigned int id : object.HAD_2_particles_ID)
      {
        HAD_2_particles_ID.push_back(id);
      }
      for(double energy : object.EM_1_energies)
      {
        EM_1_energies.push_back(energy);
      }
      for(double energy : object.EM_2_energies)
      {
        EM_2_energies.push_back(energy);
      }
      for(double energy : object.HAD_1_energies)
      {
        HAD_1_energies.push_back(energy);
      }
      for(double energy : object.HAD_2_energies)
      {
        HAD_2_energies.push_back(energy);
      }
      object.is_on = false;
      object.EM_1_particles_ID.clear();
      object.EM_2_particles_ID.clear();
      object.HAD_1_particles_ID.clear();
      object.HAD_2_particles_ID.clear();
      object.EM_1_energies.clear();
      object.EM_2_energies.clear();
      object.HAD_1_energies.clear();
      object.HAD_2_energies.clear();
      object.particles_detected.clear();
      return *this;
    }

    // Setters
    void set_is_on(bool is_on_set) {is_on = is_on_set;}
    // Getters
    vector<unsigned int> get_EM1_IDs() const {return EM_1_particles_ID;}
    vector<unsigned int> get_EM2_IDs() const {return EM_2_particles_ID;}
    vector<unsigned int> get_HAD1_IDs() const {return HAD_1_particles_ID;}
    vector<unsigned int> get_HAD2_IDs() const {return HAD_2_particles_ID;}
    vector<double> get_EM_1_energies() const {return EM_1_energies;}
    vector<double> get_EM_2_energies() const {return EM_2_energies;}
    vector<double> get_HAD_1_energies() const {return HAD_1_energies;}
    vector<double> get_HAD_2_energies() const {return HAD_2_energies;}
    vector<std::shared_ptr<Particle>> get_particles_detected() {return particles_detected;}

    // Make deafult
    void clear_calorimeter()
    {
      EM_1_particles_ID.clear();
      EM_2_particles_ID.clear();
      HAD_1_particles_ID.clear();
      HAD_1_particles_ID.clear();
      EM_1_energies.clear();
      EM_2_energies.clear();
      HAD_1_energies.clear();
      HAD_2_energies.clear();
      particles_detected.clear();
    }
    // Print Calorimeter
    void print_name()
    {
      cout<<"Calorimeter: [total particles detected, EM_1, EM_2, HAD_1, HAD_2] = ["<<EM_1_particles_ID.size()<<", "<<EM_1_particles_ID.size()<<", "<<EM_2_particles_ID.size()<<", "<<HAD_1_particles_ID.size()<<", "<<HAD_2_particles_ID.size()<<"]"<<endl;
    }
    //Detection Process


    bool detect(std::shared_ptr<Particle> & particle) // WILL RETURN TRUE FOR FINISHED AND FALSE FOR NOTHING DETECTED
    {
      if(calorimeter_efficiency  < particle->get_calorimeter_seed()) // This simulates bad behaviour of the detector
      {
        particle->change_seed_calorimeter();
        return false;
      }
      if(particle->get_mass() == 105.66) // Muons do not leave energy in calorimeter
      {
        particle->change_seed_calorimeter();
        return false;
      }
      if(particle->get_mass() == 0 and particle->get_spin() == 0.5) // Neutrinos are transparent to most detectors
      {
        particle->change_seed_calorimeter();
        return false;
      }

      // Now that the filters have finished we can detect: 
      particles_detected.push_back(particle);

      // First EM layer: 
      EM_1_particles_ID.push_back(particle->get_identifier());
      EM_1_energies.push_back(particle->energy_layer_n(0));

      // Second EM layer:
      if((particle->get_charge() == 0 and particle->get_spin() == 1 and particle->get_mass() == 0) or (particle->get_mass() == 0.511 and (particle->get_charge() == 1 or particle->get_charge()==-1)))
      {
        EM_2_particles_ID.push_back(particle->get_identifier());
        EM_2_energies.push_back(particle->energy_layer_n(1));
        particle->change_seed_calorimeter();
        return true;
      } else
      {
        EM_2_particles_ID.push_back(particle->get_identifier());
        EM_2_energies.push_back(particle->energy_layer_n(1));
      }

      // First HAD layer: (We should't have photons or electrons at this point)
      HAD_1_particles_ID.push_back(particle->get_identifier());
      HAD_1_energies.push_back(particle->energy_layer_n(2));

      // Second HAD layer:
      HAD_2_particles_ID.push_back(particle->get_identifier());
      HAD_2_energies.push_back(particle->energy_layer_n(3));

      particle->change_seed_calorimeter();
      return true;
    }
};

class MuonChamber
{
  protected:
    bool is_on{true};
    vector<std::shared_ptr<Particle>> INNER_MUON_LAYER_particles; // If particle is contained on either of this then we consider it a detected particle.
    vector<std::shared_ptr<Particle>> OUTER_MUON_LAYER_particles;
    vector<double> INNER_MUON_LAYER_energies; // probablemente innecesario
    vector<double> OUTER_MUON_LAYER_energies;
    //unsigned int number_of_muons{0};
    vector<std::shared_ptr<Particle>> particles_detected;
  
  public:
    MuonChamber() = default;
    MuonChamber(bool is_on_initial) {is_on = is_on_initial;}
    ~MuonChamber() {cout<<"Muon Chamber Destroyed"<<endl;}
    // Copy Constructor
    MuonChamber(MuonChamber const & object)
    {
      is_on = object.is_on;
      INNER_MUON_LAYER_particles.clear();
      INNER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_particles.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.INNER_MUON_LAYER_particles)
      {
        INNER_MUON_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_MUON_LAYER_particles)
      {
        OUTER_MUON_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_MUON_LAYER_energies)
      {
        INNER_MUON_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_MUON_LAYER_energies)
      {
        OUTER_MUON_LAYER_energies.push_back(energy);
      }
    }
    // Move Constructor
    MuonChamber(MuonChamber && object)
    {
      is_on = object.is_on;
      INNER_MUON_LAYER_particles.clear();
      INNER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_particles.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.INNER_MUON_LAYER_particles)
      {
        INNER_MUON_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_MUON_LAYER_particles)
      {
        OUTER_MUON_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_MUON_LAYER_energies)
      {
        INNER_MUON_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_MUON_LAYER_energies)
      {
        OUTER_MUON_LAYER_energies.push_back(energy);
      }
      // Now we free
      object.is_on = false;
      object.INNER_MUON_LAYER_particles.clear();
      object.INNER_MUON_LAYER_energies.clear();
      object.OUTER_MUON_LAYER_energies.clear();
      object.OUTER_MUON_LAYER_particles.clear();
      object.particles_detected.clear();
    }
    // Assignment operator
    MuonChamber & operator=(MuonChamber const & object)
    {
      if(&object == this) return *this;
      is_on = object.is_on;
      INNER_MUON_LAYER_particles.clear();
      INNER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_particles.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.INNER_MUON_LAYER_particles)
      {
        INNER_MUON_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_MUON_LAYER_particles)
      {
        OUTER_MUON_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_MUON_LAYER_energies)
      {
        INNER_MUON_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_MUON_LAYER_energies)
      {
        OUTER_MUON_LAYER_energies.push_back(energy);
      }
      return *this;
    }
    // Move Assignment Operator
    MuonChamber & operator=(MuonChamber && object)
    {
      std::swap(is_on, object.is_on);
      INNER_MUON_LAYER_particles.clear();
      INNER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_particles.clear();
      particles_detected.clear();
      for(std::shared_ptr<Particle> particle : object.particles_detected)
      {
        particles_detected.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.INNER_MUON_LAYER_particles)
      {
        INNER_MUON_LAYER_particles.push_back(particle);
      }
      for(std::shared_ptr<Particle> particle : object.OUTER_MUON_LAYER_particles)
      {
        OUTER_MUON_LAYER_particles.push_back(particle);
      }
      for(double energy : object.INNER_MUON_LAYER_energies)
      {
        INNER_MUON_LAYER_energies.push_back(energy);
      }
      for(double energy : object.OUTER_MUON_LAYER_energies)
      {
        OUTER_MUON_LAYER_energies.push_back(energy);
      }
      object.INNER_MUON_LAYER_particles.clear();
      object.INNER_MUON_LAYER_energies.clear();
      object.OUTER_MUON_LAYER_energies.clear();
      object.OUTER_MUON_LAYER_particles.clear();
      object.particles_detected.clear();
      return *this;
    }

    // Setters
    void set_is_on(bool is_on_set) {is_on = is_on_set;}
    // Getters
    vector<std::shared_ptr<Particle>> get_particles_detected() const {return particles_detected;}
    vector<std::shared_ptr<Particle>> get_INNER_MUON_LAYER_particles() const { return INNER_MUON_LAYER_particles;}
    vector<std::shared_ptr<Particle>> get_OUTER_MUON_LAYER_particles() const {return OUTER_MUON_LAYER_particles;}
    vector<double> get_INNER_MUON_LAYER_energies() const {return INNER_MUON_LAYER_energies;}
    vector<double> get_OUTER_MUON_LAYER_energies() const {return OUTER_MUON_LAYER_energies;}

    // Print MuonChamber
    void print_name()
    {
      cout<<"Muon Chamber: [total particles detected, Inner Layer, Outer Layer] = ["<<particles_detected.size()<<", "<<INNER_MUON_LAYER_particles.size()<<", "<<OUTER_MUON_LAYER_particles.size()<<"]"<<endl;
    }

    // Detection process
    bool detect(std::shared_ptr<Particle> & particle)
    {
      if(particle->get_mass() == 105.66) // Only Muons are detected
      {
        
        int sum_detection{0};
        // Inner layer:
        if(particle->get_muon_seeds()[0] > 10 - muon_chamber_efficiency)
        {
          INNER_MUON_LAYER_particles.push_back(particle);
          INNER_MUON_LAYER_energies.push_back(particle->energy_layer_n(7));
          sum_detection = sum_detection + 1;
        }

        // Outer layer:
        if(particle->get_muon_seeds()[1] > 10 - muon_chamber_efficiency)
        {
          OUTER_MUON_LAYER_particles.push_back(particle);
          OUTER_MUON_LAYER_energies.push_back(particle->energy_layer_n(8));
          sum_detection = sum_detection + 1;
        }

        particle->change_seeds_muon_chamber();
        // Corroborate detection
        if(sum_detection > 0) 
        {
          particles_detected.push_back(particle);
          return true;
        } else
        {
          return false;
        }
      } else
      {
        particle->change_seeds_muon_chamber();
        return false;
      }
    }

    void clear_muon_chamber()
    {
      INNER_MUON_LAYER_energies.clear();
      OUTER_MUON_LAYER_energies.clear();
      INNER_MUON_LAYER_particles.clear();
      OUTER_MUON_LAYER_particles.clear();
      particles_detected.clear();
    }
};

// Calculating Correct Detections over Total detections for each subdetector... (TRUE EFFICIENCY)
vector<double> correct_detections_over_total_detections()
{
  vector<double> results;
  int id_check_point{ids};
  // TRACKER
  Tracker MC1;
  vector<std::shared_ptr<Particle>> evento;
  vector<std::shared_ptr<Particle>> evento_calorimeter;
  vector<std::shared_ptr<Particle>> evento_muon_chamber;
  int tries{10000};
  int mal_detectados{0};
  int bien_detectados{0};
  int detections{0};

  for(int i{0}; i < tries; i++)
  {
    if(i < (int)tries/2)
    {
      std::shared_ptr<Neutron> new_particle = std::make_shared<Neutron>();
      evento.push_back(new_particle);
    } else
    {
      std::shared_ptr<Electron> new_particle = std::make_shared<Electron>();
      evento.push_back(new_particle);
    }
    std::shared_ptr<Electron> new_particle2 = std::make_shared<Electron>();
    std::shared_ptr<Muon> new_particle3 = std::make_shared<Muon>();
    evento_calorimeter.push_back(new_particle2);
    evento_muon_chamber.push_back(new_particle3);
  }

  for(int i{0}; i < tries; i++)
  {
    if(MC1.detect(evento[i]) == true)
    {
      detections = detections + 1;
      if(evento[i]->get_charge() == 0)
      {
        mal_detectados = mal_detectados + 1;
      } else
      {
        bien_detectados = bien_detectados + 1;
      }
    }
  }
  results.push_back((float)bien_detectados / (float)detections);
  
  // CALORIMETER
  detections = 0;
  Calorimeter C1;
  for(int i{0}; i < tries; i++)
  {
    if(C1.detect(evento_calorimeter[i]) == true)
    {
      detections = detections + 1;
    }
  }
  results.push_back((float)detections / (float)tries);

  //MUON CHAMBER
  detections = 0;
  MuonChamber C;

  for(int i{0}; i < tries; i++)
  {
    if(C.detect(evento_muon_chamber[i]) == true)
    {
      detections = detections + 1;
    }
  }
  results.push_back((float)detections/(float)tries);
  cout<<endl;
  ids = id_check_point;
  return results;
}

class Detector
{
  protected:
    bool is_on {true};
    int total_particles_detected{0}; // this resets
    bool has_detected {false}; // this resets
    int number_of_particles_in_event{0}; // this resets
    // Does not reset
    int num_of_incorrect_muons{0};
    int num_of_incorrect_protons{0};
    int num_of_incorrect_photons{0};
    int num_of_incorrect_neutrons{0};
    int num_of_incorrect_electrons{0};
    int num_of_correct_muons{0};
    int num_of_correct_electrons{0};
    int num_of_correct_protons{0};
    int num_of_correct_neutrons{0};
    int num_of_correct_photons{0};
    int num_of_neutrinos{0};
    int num_of_ghosts{0};

    // Energies
    vector<double> energies_of_all_particles_detected; // This resets every event as most vectors
    vector<double> energies_from_tracker; // this resets AUN NO
    vector<double> energies_from_muon_chamber; // this resets AUN NO
    vector<double> energies_from_calorimeter; // this resets AUN NO
    vector<double> hidden_energy; // this resets

    vector<std::shared_ptr<Particle>> all_particles_detected; // this resets
    vector<std::shared_ptr<Particle>> all_correct_particles; // this resets
    vector<std::shared_ptr<Particle>> all_incorrect_particles; // this resets

    // From now on we add vectors of the subdetectors.
    vector<std::shared_ptr<Tracker>> trackers_in_detector; // each subdetector calls a reset
    vector<std::shared_ptr<Calorimeter>> calorimeters_in_detector; // each subdetector calls a reset
    vector<std::shared_ptr<MuonChamber>> muon_chambers_in_detector; // each subdetector calls a reset

    // Correct Detections
    vector<std::shared_ptr<Particle>> correct_muons; // this resets
    vector<std::shared_ptr<Particle>> correct_protons; // this resets
    vector<std::shared_ptr<Particle>> correct_neutrons; // this resets
    vector<std::shared_ptr<Particle>> correct_photons; // this resets
    vector<std::shared_ptr<Particle>> correct_electrons; // this resets

    // Incorrect detection
    vector<std::shared_ptr<Particle>> incorrect_muons; // this resets
    vector<std::shared_ptr<Particle>> incorrect_protons; // this resets
    vector<std::shared_ptr<Particle>> incorrect_neutrons; // this resets
    vector<std::shared_ptr<Particle>> incorrect_photons; // this resets
    vector<std::shared_ptr<Particle>> incorrect_electrons; // this resets

    // Ghost Particles and neutrinos
    vector<std::shared_ptr<Particle>> ghost_particles; // this resets
    vector<std::shared_ptr<Particle>> neutrinos_through_detector; // this resets

    // Not Reseting Vectors
    vector<vector<std::shared_ptr<Particle>>> particles_per_event; // (Detected)
    vector<vector<std::shared_ptr<Particle>>> correct_particles_per_event;
    vector<vector<std::shared_ptr<Particle>>> incorrect_particles_per_event;
    vector<vector<std::shared_ptr<Particle>>> neutrinos_per_event;
    vector<vector<std::shared_ptr<Particle>>> ghosts_per_event;
    vector<vector<double>> energies_of_all_per_event;
    vector<vector<double>> hidden_energy_per_event;
    vector<vector<double>> energies_from_tracker_per_event;
    vector<vector<double>> energies_from_calorimeter_per_event;
    vector<vector<double>> energies_from_muon_chamber_per_event;

  public:
    Detector()
    {
      std::shared_ptr<Tracker> tracker1_ptr = std::make_shared<Tracker>(1);
      std::shared_ptr<Calorimeter> calorimeter1_ptr =  std::make_shared<Calorimeter>(1);
      std::shared_ptr<MuonChamber> muon_chamber1_ptr =  std::make_shared<MuonChamber>(1);
      trackers_in_detector.push_back(tracker1_ptr);
      calorimeters_in_detector.push_back(calorimeter1_ptr);
      muon_chambers_in_detector.push_back(muon_chamber1_ptr);
    }
    Detector(int total_particles_detected_initial) // LUEGO VER LA POSIBILIDAD DE HACER PARAMETRIZADA LA CANTIDAD DE SUBDETETCTORES.
    {
      if(total_particles_detected_initial >=0)
      {
        total_particles_detected = total_particles_detected_initial;
      } else
      {
        cout<<"Total number of particles detected must be positive. Taking absolute value."<<endl;
        total_particles_detected = total_particles_detected_initial*(-1);
      }
      std::shared_ptr<Tracker> tracker1_ptr = std::make_shared<Tracker>(1);
      std::shared_ptr<Calorimeter> calorimeter1_ptr =  std::make_shared<Calorimeter>(1);
      std::shared_ptr<MuonChamber> muon_chamber1_ptr =  std::make_shared<MuonChamber>(1);
      trackers_in_detector.push_back(tracker1_ptr);
      calorimeters_in_detector.push_back(calorimeter1_ptr);
      muon_chambers_in_detector.push_back(muon_chamber1_ptr);
    }
    ~Detector() {cout<<"Detector destroyed"<<endl;}

    // Setters
    void set_total_particles_detected(const int total_particles_detected_set) 
    {
      if(total_particles_detected_set >= 0)
      {
        total_particles_detected = total_particles_detected_set;
      } else
      {
        cout<<"Total number of particles detected must be positive. Taking absolute value."<<endl;
        total_particles_detected = total_particles_detected_set;
      }
    }
    void set_is_on(const bool is_on_set) {is_on = is_on_set;}
    void add_trackers(int number_to_add)
    {
      for(int n{0}; n < number_to_add; n++)
      {
        std::shared_ptr<Tracker> new_tracker =  std::make_shared<Tracker>();
        trackers_in_detector.push_back(new_tracker);
      }
    }
    void add_calorimeters(int number_to_add)
    {
      for(int n{0}; n < number_to_add; n++)
      {
        std::shared_ptr<Calorimeter> new_calorimeter =  std::make_shared<Calorimeter>();
        calorimeters_in_detector.push_back(new_calorimeter);
      }
    }
    void add_muon_chambers(int number_to_add)
    {
      for(int n{0}; n < number_to_add; n++)
      {
        std::shared_ptr<MuonChamber> new_MC =  std::make_shared<MuonChamber>();
        muon_chambers_in_detector.push_back(new_MC);
      }
    }
    
    // Getters
    int get_total_particles_detected() const {return total_particles_detected;}
    vector<std::shared_ptr<Tracker>> get_trackers_in_detector() const {return trackers_in_detector;}
    vector<std::shared_ptr<Calorimeter>> get_calorimeters_in_detector() const {return calorimeters_in_detector;}
    vector<std::shared_ptr<MuonChamber>> get_muon_chambers_in_detector() const {return muon_chambers_in_detector;}
    vector<vector<std::shared_ptr<Particle>>> get_correct_detections() const
    {
      vector<vector<std::shared_ptr<Particle>>> vector_of_vectors;
      vector_of_vectors.push_back(correct_electrons);
      vector_of_vectors.push_back(correct_muons);
      vector_of_vectors.push_back(correct_photons);
      vector_of_vectors.push_back(correct_protons);
      vector_of_vectors.push_back(correct_neutrons);
      return  vector_of_vectors;
    }
    vector<vector<std::shared_ptr<Particle>>> get_all_particles_per_event() const {return particles_per_event;}
    vector<vector<std::shared_ptr<Particle>>> get_correct_particles_per_event() const {return correct_particles_per_event;}
    vector<vector<std::shared_ptr<Particle>>> get_incorrect_particles_per_event() const {return incorrect_particles_per_event;}
    
    // Detector Print
    void print_name()
    {
      cout<<"... PRINTING DETECTOR INFORMATION ..."<<endl;
      cout<<"Detector: [Trackers, Calorimeters, MuonChambers] = ["<<trackers_in_detector.size()<<", "<<calorimeters_in_detector.size()<<", "<<muon_chambers_in_detector.size()<<"]"<<endl;
      cout<<"Number of Particles detected: "<<total_particles_detected<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"Protons Detected: "<<num_of_correct_protons + num_of_incorrect_protons<<" Correct detections: "<<num_of_correct_protons<<endl;
      cout<<"Electrons Detected: "<<num_of_correct_electrons + num_of_incorrect_electrons<<" Correct detections: "<<num_of_correct_electrons<<endl;
      cout<<"Muons Detected: "<<num_of_correct_muons + num_of_incorrect_muons<<" Correct detections: "<<num_of_correct_muons<<endl;
      cout<<"Photons Detected: "<<num_of_correct_photons + num_of_incorrect_photons<<" Correct detections: "<<num_of_correct_photons<<endl;
      cout<<"Neutrons Detected: "<<num_of_correct_neutrons + num_of_incorrect_neutrons<<" Correct detections: "<<num_of_correct_neutrons<<endl;
      cout<<"Number of not recognized particles: "<<num_of_ghosts<<endl;
    }

    // Detection Process
    // Here we will iterate on each list of trackers, muon chambers and calorimeters and we will perform a detection
    // dependindg on the nature of the subdetetctor we can actually call functions of detection from the specific class.
    bool detect(vector<std::shared_ptr<Particle>> event)
    {    
      if(is_on == false)
      {
        cout<<" ----> Detector is off."<<endl;
        return false;
      }

      // Calculate true efficiencies of the detector:
      vector<double> true_efficiencies{correct_detections_over_total_detections()};
      cout<<endl;
      cout<<"...Detection Started..."<<endl;
      cout<<"------------------------"<<endl;
      double true_efficiency_tracker{true_efficiencies[0]};
      double true_efficiency_calorimeter{true_efficiencies[1]};
      double true_efficiency_muon_chamber{true_efficiencies[2]};
      cout<<"True Efficiency of Tracker: "<<true_efficiency_tracker<<endl;
      cout<<"True Efficiency of Calorimeter: "<<true_efficiency_calorimeter<<endl;
      cout<<"True Efficiency of Muon Chamber: "<<true_efficiency_muon_chamber<<endl; cout<<endl;
      //

      number_of_particles_in_event = event.size();
      // Trackers
      for(int i{0}; i < trackers_in_detector.size(); i++) // Repeat for each Tracker
      {
        cout<<"Tracker Number: "<<i+1<<endl;
        cout<<"--------------"<<endl;
        for (int j{0}; j < event.size(); j++) // Repeat for each particle in the event
        {
          cout<<"Particle ID: "<<event[j]->get_identifier()<<endl;
          cout<<"---"<<endl;
          // Need to treat Tau separately
          if(event[j]->get_mass() == 1776.86)
          {
            if(event[j]->get_decayed_yet() == false) // Extra security, since it is already set in the decay function
            { 
              event[j]->decay();
            } // if already decayed, then it does not decay again
            if(event[j]->get_decay_type() == "leptonic")
            {
              number_of_particles_in_event = number_of_particles_in_event + 2;
              trackers_in_detector[i]->detect(event[j]->get_decay_products()[0]);
              trackers_in_detector[i]->detect(event[j]->get_decay_products()[1]);
              trackers_in_detector[i]->detect(event[j]->get_decay_products()[2]);
            } else // if decay is hadronic
            {
              trackers_in_detector[i]->detect(event[j]); // Will interact as proton.
            }
          } else // if not a tau particle
          {
            trackers_in_detector[i]->detect(event[j]);
          }
          cout<<"Finished with Particle: "<<event[j]->get_identifier()<<endl;
          cout<<"-"<<endl;
        }
        // Print result for Tracker
        trackers_in_detector[i]->print_name();
      } // Trackers finished.
      cout<<"Finished with Trackers..."<<endl; cout<<endl;

      //Calorimeters
      for(int i{0}; i < calorimeters_in_detector.size(); i++) // Repeat for each Calorimeter
      {
        cout<<"Calorimeter Number: "<<i+1<<endl;
        cout<<"--------------"<<endl;
        for(int j{0}; j < event.size(); j++) // Repeat for each particle in the event
        {
          cout<<"Particle ID: "<<event[j]->get_identifier()<<endl;
          cout<<"---"<<endl;
          // Tau Treatment
          if(event[j]->get_mass() == 1776.86)
          {
            if(event[j]->get_decayed_yet() == false) // Extra security, since it is already set in the decay function
            {
              event[j]->decay();
            } // if already decayed, then it does not decay again
            if(event[j]->get_decay_type() == "leptonic")
            {
              calorimeters_in_detector[i]->detect(event[j]->get_decay_products()[0]);
              calorimeters_in_detector[i]->detect(event[j]->get_decay_products()[1]);
              calorimeters_in_detector[i]->detect(event[j]->get_decay_products()[2]);
            } else // if tau decays hadronically
            {
              calorimeters_in_detector[i]->detect(event[j]);
            }
          } else // if not a Tau
          {
            calorimeters_in_detector[i]->detect(event[j]);
          }
          cout<<"Finished with Particle: "<<event[j]->get_identifier()<<endl;
          cout<<"-"<<endl;
        }
        // Print result for Calorimeter
        calorimeters_in_detector[i]->print_name();
      } // Calorimeters finished.
      cout<<"Finished with calorimeters..."<<endl; cout<<endl;

      // Muon Chambers
      for(int i{0}; i < muon_chambers_in_detector.size(); i++) // Repeat for each MuonChamber
      {
        cout<<"Muon Chamber Number: "<<i+1<<endl;
        cout<<"--------------"<<endl;
        for(int j{0}; j < event.size(); j++) // For each Particle
        {
          cout<<"Particle ID: "<<event[j]->get_identifier()<<endl;
          cout<<"---"<<endl;
          // Tau Treatment
          if(event[j]->get_mass() == 1776.86)
          {
            if(event[j]->get_decayed_yet() == false) // Extra security.
            {
              event[j]->decay();
            } // if already decayed, then it does not decay again
            if(event[j]->get_decay_type() == "leptonic")
            {
              muon_chambers_in_detector[i]->detect(event[j]->get_decay_products()[0]);
              muon_chambers_in_detector[i]->detect(event[j]->get_decay_products()[1]);
              muon_chambers_in_detector[i]->detect(event[j]->get_decay_products()[2]);
            } else // if tau decays hadronically
            {
              muon_chambers_in_detector[i]->detect(event[j]);
            }
          } else // if not a tau
          {
            muon_chambers_in_detector[i]->detect(event[j]);
          }
          cout<<"Finished with Particle: "<<event[j]->get_identifier()<<endl;
          cout<<"-"<<endl;
        }
        // Print result for Muon Chamber
        muon_chambers_in_detector[i]->print_name();
      }
      cout<<endl; cout<<"...DETECTION PROCESS FINISHED..."<<endl; cout<<endl;
      has_detected = true;

      // ANALYSIS PROCESS STARTS HERE
      // ----------------------------
      // Check if the "geometry" of the detector is suitable to analyse the data
      if(trackers_in_detector.size() == 0 or calorimeters_in_detector.size() == 0 or muon_chambers_in_detector.size() == 0)
      {
        cout<<" ------> GEOMETRY OF DETECTOR NOT SUITABLE FOR PROPER ANALYSIS <------ "<<endl;
        return false;
      }
      
      // Check if the detector has detected an event
      if(has_detected == false)
      {
        cout<<"No detection yet"<<endl;
        return false;
      }
      // End of checking if there has been a detected event.

      static int events_analysed{1}; // ESTATIUCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

      cout<<"Starting Analysis of Last Event"<<endl;
      cout<<"-------------------------------"<<endl;
      for(std::shared_ptr<Particle> particle : event) // For each particle in the event
      {
        // We distinguish between Leptonic Tau and not leptonic Tau.
        if(particle->get_mass() == 1776.86)
        {
          if(particle->get_decay_type() == "leptonic")
          {
            for(int decay_index{0}; decay_index < 3; decay_index++) // For each deacy product we must do the same analysis...
            {
              // Define temporal bool variables for the particle.
              bool detected_tracker{false};
              bool detected_calorimeter{false};
              bool detected_muon_chamber{false};
              bool detected_EM{false};
              bool detected_HAD{false};
              bool detected_INNER_TRACKER{false};
              bool detected_OUTER_TRACKER{false};
              bool detected_STRIP_LAYER_TRACKER{false};
              bool detected_INNER_MUON{false};
              bool detected_OUTER_MUON{false};

              // Tracker Analysis
              //------------------
              int number_of_trackers_that_detected{0}; // temporal variable.
              int num_of_INNER_TRACKER_that_detected{0};
              int num_of_OUTER_TRACKER_that_detected{0};
              int num_of_STRIP_TRACKER_that_detected{0};
              for(std::shared_ptr<Tracker> tracker : trackers_in_detector) // For each particle we are iterating through the trackers.
              {
                for(std::shared_ptr<Particle> detected_particle : tracker->get_particles_detected()) // Iterate on list of particles detected in each of the trackers
                {
                  if(detected_particle->get_identifier() == particle->get_decay_products()[decay_index]->get_identifier())
                 {
                   number_of_trackers_that_detected = number_of_trackers_that_detected + 1;
                  }
                }
                for(std::shared_ptr<Particle> p_in_inner : tracker->get_particles_INNER_PIXEL_LAYER())
                {
                  if(p_in_inner->get_identifier() == particle->get_decay_products()[decay_index]->get_identifier())
                  {
                    num_of_INNER_TRACKER_that_detected++;
                  }
                }
                for(std::shared_ptr<Particle> p_in_outer : tracker->get_particles_OUTER_PIXEL_LAYER())
                {
                  if(p_in_outer->get_identifier() == particle->get_decay_products()[decay_index]->get_identifier())
                  {
                    num_of_OUTER_TRACKER_that_detected++;
                  }
                }
                for(std::shared_ptr<Particle> p_in_strip : tracker->get_particles_STRIP_LAYER())
                {
                  if(p_in_strip->get_identifier() == particle->get_decay_products()[decay_index]->get_identifier())
                  {
                    num_of_STRIP_TRACKER_that_detected++;
                  }
                }
              } // Finished iterating particle in each detector
              //Now we make a statistical analysis considering the efficiency of our trackers...
              if(number_of_detections_for_confidence(true_efficiency_tracker, trackers_in_detector.size()) < number_of_trackers_that_detected)
              {
                detected_tracker = true;
              } else
              {
                if(number_of_detections_for_confidence(true_efficiency_tracker, trackers_in_detector.size()) == 1 and number_of_trackers_that_detected > 0)
                {
                  detected_tracker = true;
                }
              }

              // EACH LAYER OF TRACKER TREATMENT WILL BE JUST AVERAGE DETECTIONS:
              if(num_of_INNER_TRACKER_that_detected / trackers_in_detector.size() >= 0.5 and detected_tracker == true)
              {
                detected_INNER_TRACKER = true;
              }
              if(num_of_OUTER_TRACKER_that_detected / trackers_in_detector.size() >= 0.5 and detected_tracker == true)
              {
                detected_OUTER_TRACKER = true;
              }
              if(num_of_STRIP_TRACKER_that_detected / trackers_in_detector.size() >= 0.5 and detected_tracker == true)
              {
                detected_STRIP_LAYER_TRACKER = true;
              }

              // Calorimeter Analysis
              // --------------------
              int number_of_calorimeters_that_detected{0};
              int number_of_EM_that_detected{0};
              int number_of_HAD_that_detected{0};

              for(std::shared_ptr<Calorimeter> calorimeter : calorimeters_in_detector) // For each particle we are iterating through the calorimeters.
              {
                for(std::shared_ptr<Particle> detected_particle : calorimeter->get_particles_detected()) // Iterate on list of particles detected in each of the calorimeter
                {
                  if(detected_particle->get_identifier() == particle->get_decay_products()[decay_index]->get_identifier())
                  {
                    number_of_calorimeters_that_detected = number_of_calorimeters_that_detected + 1;
                  }
                }
                for(unsigned int id : calorimeter->get_EM1_IDs()) // Iterate vector containing the id of detected particles
                {
                  if(id == particle->get_decay_products()[decay_index]->get_identifier())
                  {
                    number_of_EM_that_detected = number_of_EM_that_detected + 1;
                  }
                }
                for(unsigned int id : calorimeter->get_HAD1_IDs())
                {
                  if(id == particle->get_decay_products()[decay_index]->get_identifier())
                  {
                    number_of_HAD_that_detected = number_of_HAD_that_detected + 1;
                  }
                }
              } // Finished iterating particle in each calorimeter
              //Now we make a statistical analysis considering the efficiency of our calorimeters...
              if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) < number_of_calorimeters_that_detected)
              {
                detected_calorimeter = true;
              } else
              {
                if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) == 1 and number_of_calorimeters_that_detected > 0)
                {
                  detected_calorimeter = true;
                }
              }
              if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) < number_of_EM_that_detected)
              {
                detected_EM = true;
              } else
              {
                if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) == 1 and number_of_EM_that_detected > 0)
                {
                  detected_EM = true;
                }
              }
              if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) < number_of_HAD_that_detected)
              {
                detected_HAD = true;
              } else
              {
                if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) == 1 and number_of_HAD_that_detected > 0)
                {
                  detected_HAD = true;
                }
              }
              
              // MuonChamber Analysis
              // --------------------
              int number_of_muon_chambers_that_detected{0};
              int num_inner_muon_that_detected{0};
              int num_outer_muon_that_detected{0};
              for(std::shared_ptr<MuonChamber> muon_chamber : muon_chambers_in_detector)// For each particle we are iterating through the MuonChambers.
              {
                for(std::shared_ptr<Particle> detected_particle : muon_chamber->get_particles_detected()) // for each particle of the vector of detected particles of each muon chamber
                {
                  if(particle->get_decay_products()[decay_index]->get_identifier() == detected_particle->get_identifier())
                  {
                    number_of_muon_chambers_that_detected = number_of_muon_chambers_that_detected + 1;
                  }
                }
                for(std::shared_ptr<Particle> p_in_inner : muon_chamber->get_INNER_MUON_LAYER_particles())
                {
                  if(particle->get_decay_products()[decay_index]->get_identifier() == p_in_inner->get_identifier())
                  {
                    num_inner_muon_that_detected++;
                  }
                }
                for(std::shared_ptr<Particle> p_in_outer : muon_chamber->get_OUTER_MUON_LAYER_particles())
                {
                  if(particle->get_decay_products()[decay_index]->get_identifier() == p_in_outer->get_identifier())
                  { 
                    num_outer_muon_that_detected++;
                  }
                }
              }// Finished iterating particle in each Muon Chamber
              //Now we make a statistical analysis considering the efficiency of our MuonChambers...
              if(number_of_detections_for_confidence(true_efficiency_muon_chamber, muon_chambers_in_detector.size()) < number_of_muon_chambers_that_detected)
              {
                detected_muon_chamber = true;
              } else
              {
                if(number_of_detections_for_confidence(true_efficiency_muon_chamber, muon_chambers_in_detector.size()) == 1 and number_of_muon_chambers_that_detected > 0)
                {
                  detected_muon_chamber = true;
                }
              }

              // EACH LAYER OF MUON CHAMBER TREATMENT WILL BE JUST AVERAGE DETECTIONS:
              if(num_inner_muon_that_detected / muon_chambers_in_detector.size() >= 0.5 and detected_muon_chamber == true)
              {
                detected_INNER_MUON = true;
              }
              if(num_outer_muon_that_detected / muon_chambers_in_detector.size() >= 0.5 and detected_muon_chamber == true)
              {
                detected_OUTER_MUON = true;
              }
              
              // Maybe add some error in identification--------
              if(rand() % 10 + 1 > detector_analysis_efficiency)
              {
                cout<<"Random Error identifying activated"<<endl;
                detected_tracker = rand() % 2 == 0;
                detected_EM = rand() % 2 == 0;
                detected_HAD = rand() % 2 == 0;
                detected_INNER_MUON = rand() % 2 == 0;
                detected_muon_chamber = rand() % 2 == 0;
                detected_INNER_TRACKER = rand() % 2 == 0;
                detected_OUTER_MUON = rand() % 2 == 0;
                detected_OUTER_TRACKER= rand() % 2 == 0;
                detected_STRIP_LAYER_TRACKER = rand() % 2 == 0;

                if(detected_EM == false and detected_HAD == false)
                {
                  detected_calorimeter = false;
                }
                if(detected_INNER_MUON == false and detected_OUTER_MUON == false) 
                {
                  detected_muon_chamber = false;
                }
                if(detected_INNER_TRACKER == false and detected_OUTER_TRACKER == false and detected_STRIP_LAYER_TRACKER == false)
                {
                  detected_tracker = false;
                }
              }
              // ----------------------------------------------

              // AT THIS POINT WE HAVE ALL THE BOOLS WE NEED TO START ANALYSING THE PARTICLE
              cout<<"Particle ID: "<<particle->get_decay_products()[decay_index]->get_identifier()<<endl;
              cout<<"------------"<<endl;
              cout<<"Detections with at least 0.95 confidence: "<<"Tracker = "<<detected_tracker<<", "<<"Calorimeter = [Overall, EM, HAD]: ["<<detected_calorimeter<<", "<<detected_EM<<", "<<detected_HAD<<"], MuonChamber =  "<<detected_muon_chamber<<endl;
              cout<<"-----------------------------------------"<<endl;
              cout<<endl;
              
              if(detected_tracker == true or detected_calorimeter == true or detected_muon_chamber == true)
              {
                total_particles_detected = total_particles_detected + 1;
                all_particles_detected.push_back(particle->get_decay_products()[decay_index]);
              }

              // NOW WE IDENTIFY THE PARTICLE
              bool identified_yet{false};
              double energy_detected{0};
              double energy_from_tracker{0};
              double energy_from_muon_chamber{0};
              double energy_from_calorimeter{0};
              // ENERGY TREATMENT DETECTION
              if(detected_EM == true)
              {
                energy_detected = energy_detected + particle->get_decay_products()[decay_index]->energy_layer_n(0) + particle->get_decay_products()[decay_index]->energy_layer_n(1);
                energy_from_calorimeter = energy_from_calorimeter + particle->get_decay_products()[decay_index]->energy_layer_n(0) + particle->get_decay_products()[decay_index]->energy_layer_n(1);
              }
              if(detected_HAD == true)
              {
                energy_detected = energy_detected + particle->get_decay_products()[decay_index]->energy_layer_n(2) + particle->get_decay_products()[decay_index]->energy_layer_n(3);
                energy_from_calorimeter = energy_from_calorimeter + particle->get_decay_products()[decay_index]->energy_layer_n(2) + particle->get_decay_products()[decay_index]->energy_layer_n(3);
              }
              if(detected_INNER_TRACKER == true)
              {
                energy_detected = energy_detected + particle->get_decay_products()[decay_index]->energy_layer_n(4);
                energy_from_tracker = energy_from_tracker + particle->get_decay_products()[decay_index]->energy_layer_n(4);
              }
              if(detected_OUTER_TRACKER == true)
              {
                energy_detected = energy_detected + particle->get_decay_products()[decay_index]->energy_layer_n(5);
                energy_from_tracker = energy_from_tracker + particle->get_decay_products()[decay_index]->energy_layer_n(5);
              }
              if(detected_STRIP_LAYER_TRACKER == true)
              {
                energy_detected = energy_detected + particle->get_decay_products()[decay_index]->energy_layer_n(6);
                energy_from_tracker = energy_from_tracker + particle->get_decay_products()[decay_index]->energy_layer_n(6);
              }
              if(detected_INNER_MUON == true)
              {
                energy_detected = energy_detected + particle->get_decay_products()[decay_index]->energy_layer_n(7);
                energy_from_muon_chamber = energy_from_muon_chamber + particle->get_decay_products()[decay_index]->energy_layer_n(7);
              }
              if(detected_OUTER_MUON)
              {
                energy_detected = energy_detected + particle->get_decay_products()[decay_index]->energy_layer_n(8);
                energy_from_muon_chamber = energy_from_muon_chamber + particle->get_decay_products()[decay_index]->energy_layer_n(8);
              }
              energies_of_all_particles_detected.push_back(energy_detected);
              energies_from_calorimeter.push_back(energy_from_calorimeter);
              energies_from_muon_chamber.push_back(energy_from_muon_chamber);
              energies_from_tracker.push_back(energy_from_tracker);
              // Finished with energy treatment

              if(detected_tracker == false and detected_EM == true and detected_muon_chamber == false and detected_HAD == false)
              {
                // We think this particle could be a Photon so we check against the real information and then add the particle to
                // list of correctly detected photons or incorrectly detected photon
                if(particle->get_decay_products()[decay_index]->get_mass() == 0 and particle->get_decay_products()[decay_index]->get_charge() == 0 and particle->get_decay_products()[decay_index]->get_spin() == 1)
                {
                  identified_yet = true;
                  correct_photons.push_back(particle->get_decay_products()[decay_index]);
                  all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                } else
                {
                  identified_yet = true;
                  incorrect_photons.push_back(particle->get_decay_products()[decay_index]);
                  all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                }
              } // Finished checking for Photon
              if(detected_tracker == true)
              {
                if(detected_EM == true and detected_HAD == false) // For Electron
                {
                  // We think this particle must be an electron
                  if(particle->get_decay_products()[decay_index]->get_mass() == 0.511)
                  {
                    identified_yet = true;
                    correct_electrons.push_back(particle->get_decay_products()[decay_index]);
                    all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                  } else
                  {
                    identified_yet = true;
                    incorrect_electrons.push_back(particle->get_decay_products()[decay_index]);
                    all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                  }
                }
                if(detected_EM == true and detected_HAD == true) // For Proton
                {
                  // We think this must be a Proton
                  if(particle->get_decay_products()[decay_index]->get_mass() == 938)
                  {
                    identified_yet = true;
                    correct_protons.push_back(particle->get_decay_products()[decay_index]);
                    all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                  } else
                  {
                    identified_yet = true;
                    incorrect_protons.push_back(particle->get_decay_products()[decay_index]);
                    all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                  }
                }
                if(detected_calorimeter == false and detected_muon_chamber == true) // For Muon
                {
                  // We think this must be a muon
                  if(particle->get_decay_products()[decay_index]->get_mass() == 105.66)
                  {
                    identified_yet = true;
                    correct_muons.push_back(particle->get_decay_products()[decay_index]);
                    all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                  } else
                  {
                    identified_yet = true;
                    incorrect_muons.push_back(particle->get_decay_products()[decay_index]);
                    all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                  }
                }

                if(detected_calorimeter == false and detected_muon_chamber == false) // Only Tracker Signature
                {
                  cout<<"........................................"<<endl;
                  cout<<"Only Tracker made a detection, particle can be either Muon, Electron or Proton"<<endl;
                  cout<<"........................................"<<endl;
                  int random{rand() % 3 + 1};
                  if(random == 1)
                  {
                    // Electron Case
                    if(particle->get_decay_products()[decay_index]->get_mass() == 0.511)
                    {
                      identified_yet = true;
                      correct_electrons.push_back(particle->get_decay_products()[decay_index]);
                      all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                    } else
                    {
                      identified_yet = true;
                      incorrect_electrons.push_back(particle->get_decay_products()[decay_index]);
                      all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                    }
                  }
                  if(random == 2)
                  {
                    // Muon case
                    if(particle->get_decay_products()[decay_index]->get_mass() == 105.66)
                    {
                      identified_yet = true;
                      correct_muons.push_back(particle->get_decay_products()[decay_index]);
                      all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                    } else
                    {
                      identified_yet = true;
                      incorrect_muons.push_back(particle->get_decay_products()[decay_index]);
                      all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                    }
                  }
                  if(random == 3)
                  {
                    // Proton case
                    if(particle->get_decay_products()[decay_index]->get_mass() == 938)
                    {
                      identified_yet = true;
                      correct_protons.push_back(particle->get_decay_products()[decay_index]);
                      all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                    } else
                    {
                      identified_yet = true;
                      incorrect_protons.push_back(particle->get_decay_products()[decay_index]);
                      all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                    }
                  }
                }
              } // Finished distinguishing tracker signatures
              if(detected_tracker == false and detected_muon_chamber == false and detected_EM == true and detected_HAD == true)
              {
                // We think this must be a neutron
                if(particle->get_decay_products()[decay_index]->get_mass() == 939.5)
                {
                  identified_yet = true;
                  correct_neutrons.push_back(particle->get_decay_products()[decay_index]);
                  all_correct_particles.push_back(particle->get_decay_products()[decay_index]);
                } else
                {
                  identified_yet = true;
                  incorrect_neutrons.push_back(particle->get_decay_products()[decay_index]);
                  all_incorrect_particles.push_back(particle->get_decay_products()[decay_index]);
                }
              } // Finished looking for Neutron
              if(identified_yet == false) // IF PARTICLE STILL NOT IDENTIFIED THEN IT GOES TO GHOST PARTICLE VECTOR
              {
                if(particle->get_decay_products()[decay_index]->get_mass() == 0 and particle->get_decay_products()[decay_index]->get_spin() == 0.5)
                {
                  neutrinos_through_detector.push_back(particle->get_decay_products()[decay_index]);
                  hidden_energy.push_back(particle->get_decay_products()[decay_index]->get_true_energy());
                } else
                {
                  ghost_particles.push_back(particle->get_decay_products()[decay_index]);
                  hidden_energy.push_back(particle->get_decay_products()[decay_index]->get_true_energy());
                }
              }

            } // FINISHED WITH DECAY PRODUCTS...
            continue;
          } 
        } // FINISHED WITH TAU
        // -------------------

        // Define temporal bool variables for the particle.
        bool detected_tracker{false};
        bool detected_calorimeter{false};
        bool detected_muon_chamber{false};
        bool detected_EM{false};
        bool detected_HAD{false};
        bool detected_INNER_TRACKER{false};
        bool detected_OUTER_TRACKER{false};
        bool detected_STRIP_LAYER_TRACKER{false};
        bool detected_INNER_MUON{false};
        bool detected_OUTER_MUON{false};

        // Tracker Analysis
        //------------------
        int number_of_trackers_that_detected{0}; // temporal variable.
        int num_of_INNER_TRACKER_that_detected{0};
        int num_of_OUTER_TRACKER_that_detected{0};
        int num_of_STRIP_TRACKER_that_detected{0};
        for(std::shared_ptr<Tracker> tracker : trackers_in_detector) // For each particle we are iterating through the trackers.
        {
          for(std::shared_ptr<Particle> detected_particle : tracker->get_particles_detected()) // Iterate on list of particles detected in each of the trackers
          {
            if(detected_particle->get_identifier() == particle->get_identifier())
            {
              number_of_trackers_that_detected = number_of_trackers_that_detected + 1;
            }
          }
          for(std::shared_ptr<Particle> p_in_inner : tracker->get_particles_INNER_PIXEL_LAYER())
          {
            if(p_in_inner->get_identifier() == particle->get_identifier())
            {
              num_of_INNER_TRACKER_that_detected++;
            }
          }
          for(std::shared_ptr<Particle> p_in_outer : tracker->get_particles_OUTER_PIXEL_LAYER())
          {
            if(p_in_outer->get_identifier() == particle->get_identifier())
            {
              num_of_OUTER_TRACKER_that_detected++;
            }
          }
          for(std::shared_ptr<Particle> p_in_strip : tracker->get_particles_STRIP_LAYER())
          {
            if(p_in_strip->get_identifier() == particle->get_identifier())
            {
              num_of_STRIP_TRACKER_that_detected++;
            }
          }
        } // Finished iterating particle in each detector
        //Now we make a statistical analysis considering the efficiency of our trackers...
        if(number_of_detections_for_confidence(true_efficiency_tracker, trackers_in_detector.size()) < number_of_trackers_that_detected)
        {
          detected_tracker = true;
        } else
        {
          if(number_of_detections_for_confidence(true_efficiency_tracker, trackers_in_detector.size()) == 1 and number_of_trackers_that_detected > 0)
          {
            detected_tracker = true;
          }
        }

        // EACH LAYER OF TRACKER TREATMENT WILL BE JUST AVERAGE DETECTIONS:
        if(num_of_INNER_TRACKER_that_detected / trackers_in_detector.size() >= 0.5 and detected_tracker == true)
        {
          detected_INNER_TRACKER = true;
        }
        if(num_of_OUTER_TRACKER_that_detected / trackers_in_detector.size() >= 0.5 and detected_tracker == true)
        {
          detected_OUTER_TRACKER = true;
        }
        if(num_of_STRIP_TRACKER_that_detected / trackers_in_detector.size() >= 0.5 and detected_tracker == true)
        {
          detected_STRIP_LAYER_TRACKER = true;
        }
        

        // Calorimeter Analysis
        // --------------------
        int number_of_calorimeters_that_detected{0};
        int number_of_EM_that_detected{0};
        int number_of_HAD_that_detected{0};
        for(std::shared_ptr<Calorimeter> calorimeter : calorimeters_in_detector) // For each particle we are iterating through the calorimeters.
        {
          for(std::shared_ptr<Particle> detected_particle : calorimeter->get_particles_detected()) // Iterate on list of particles detected in each of the calorimeter
          {
            if(detected_particle->get_identifier() == particle->get_identifier())
            {
              number_of_calorimeters_that_detected = number_of_calorimeters_that_detected + 1;
            }
          }
          for(unsigned int id : calorimeter->get_EM1_IDs()) // Iterate vector containing the id of detected particles
          {
            if(id == particle->get_identifier())
            {
              number_of_EM_that_detected = number_of_EM_that_detected + 1;
            }
          }
          for(unsigned int id : calorimeter->get_HAD1_IDs())
          {
            if(id == particle->get_identifier())
            {
              number_of_HAD_that_detected = number_of_HAD_that_detected + 1;
            }
          }
        } // Finished iterating particle in each detector
        //Now we make a statistical analysis considering the efficiency of our calorimeters...
        if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) < number_of_calorimeters_that_detected)
        {
          detected_calorimeter = true;
        } else
        {
          if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) == 1 and number_of_calorimeters_that_detected > 0)
          {
            detected_calorimeter = true;
          }
        }
        if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) < number_of_EM_that_detected)
        {
          detected_EM = true;
        } else

        {
          if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) == 1 and number_of_EM_that_detected > 0)
          {
            detected_EM = true;
          }
        }
        if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) < number_of_HAD_that_detected)
        {
          detected_HAD = true;
        } else
        {
          if(number_of_detections_for_confidence(true_efficiency_calorimeter, calorimeters_in_detector.size()) == 1 and number_of_HAD_that_detected > 0)
          {
            detected_HAD = true;
          }
        }
        

        // MuonChamber Analysis
        // --------------------
        int number_of_muon_chambers_that_detected{0};
        int num_inner_muon_that_detected{0};
        int num_outer_muon_that_detected{0};
        for(std::shared_ptr<MuonChamber> muon_chamber : muon_chambers_in_detector)// For each particle we are iterating through the MuonChambers.
        {
          for(std::shared_ptr<Particle> detected_particle : muon_chamber->get_particles_detected()) // for each particle of the vector of detected particles of each muon chamber
          {
            if(particle->get_identifier() == detected_particle->get_identifier())
            {
              number_of_muon_chambers_that_detected = number_of_muon_chambers_that_detected + 1;
            }
          }
          for(std::shared_ptr<Particle> p_in_inner : muon_chamber->get_INNER_MUON_LAYER_particles())
          {
            if(particle->get_identifier() == p_in_inner->get_identifier())
            {
              num_inner_muon_that_detected++;
            }
          }
          for(std::shared_ptr<Particle> p_in_outer : muon_chamber->get_OUTER_MUON_LAYER_particles())
          {
            if(particle->get_identifier() == p_in_outer->get_identifier())
            {
              num_outer_muon_that_detected++;
            }
          }
        }// Finished iterating particle in each detector
        //Now we make a statistical analysis considering the efficiency of our MuonChambers...
        if(number_of_detections_for_confidence(true_efficiency_muon_chamber, muon_chambers_in_detector.size()) < number_of_muon_chambers_that_detected)
        {
          detected_muon_chamber = true;
        } else
        {
          if(number_of_detections_for_confidence(true_efficiency_muon_chamber, muon_chambers_in_detector.size()) == 1 and number_of_muon_chambers_that_detected > 0)
          {
            detected_muon_chamber = true;
          }
        }

        // EACH LAYER OF MUON CHAMBER TREATMENT WILL BE JUST AVERAGE DETECTIONS:
        if(num_inner_muon_that_detected / muon_chambers_in_detector.size() >= 0.5 and detected_muon_chamber == true)
        {
          detected_INNER_MUON = true;
        }
        if(num_outer_muon_that_detected / muon_chambers_in_detector.size() >= 0.5 and detected_muon_chamber == true)
        {
          detected_OUTER_MUON = true;
        }

        // Maybe add some error in identification--------
        if(rand() % 10 + 1 > detector_analysis_efficiency)
        {
          cout<<"Random Error identifying activated"<<endl;
          detected_tracker = rand() % 2 == 0;
          detected_EM = rand() % 2 == 0;
          detected_HAD = rand() % 2 == 0;
          detected_INNER_MUON = rand() % 2 == 0;
          detected_muon_chamber = rand() % 2 == 0;
          detected_INNER_TRACKER = rand() % 2 == 0;
          detected_OUTER_MUON = rand() % 2 == 0;
          detected_OUTER_TRACKER= rand() % 2 == 0;
          detected_STRIP_LAYER_TRACKER = rand() % 2 == 0;

          if(detected_EM == false and detected_HAD == false)
          {
            detected_calorimeter = false;
          }
          if(detected_INNER_MUON == false and detected_OUTER_MUON == false) 
          {
            detected_muon_chamber = false;
          }
          if(detected_INNER_TRACKER == false and detected_OUTER_TRACKER == false and detected_STRIP_LAYER_TRACKER == false)
          {
            detected_tracker = false;
          }
        }
        // ----------------------------------------------

        // AT THIS POINT WE HAVE ALL THE BOOLS WE NEED TO START ANALYSING THE PARTICLE
        cout<<"Particle ID: "<<particle->get_identifier()<<endl;
        cout<<"------------"<<endl;
        cout<<"Detections 0.95 confidence: "<<"Tracker = [Overall, INNER, OUTER, STRIP]: ["<<detected_tracker<<", "<<detected_INNER_TRACKER<<", "<<detected_OUTER_TRACKER<<", "<<detected_STRIP_LAYER_TRACKER<<"] Calorimeter = [Overall, EM, HAD]: ["<<detected_calorimeter<<", "<<detected_EM<<", "<<detected_HAD<<"], MuonChamber =  [Overall, INNER, OUTER]: ["<<detected_muon_chamber<<", "<<detected_INNER_MUON<<", "<<detected_OUTER_MUON<<"]"<<endl;
        cout<<"-----------------------------------------"<<endl;
        cout<<endl;

        if(detected_tracker == true or detected_calorimeter == true or detected_muon_chamber == true)
        {
          total_particles_detected = total_particles_detected + 1;
          all_particles_detected.push_back(particle); // I add the particle to the detected particles of this event
        }

        // NOW WE IDENTIFY THE PARTICLE
        bool identified_yet{false};
        
        double energy_detected{0};
        double energy_from_tracker{0};
        double energy_from_muon_chamber{0};
        double energy_from_calorimeter{0};
        // ENERGY TREATMENT DETECTION
        if(detected_EM == true)
        {
          energy_detected = energy_detected + particle->energy_layer_n(0) + particle->energy_layer_n(1);
          energy_from_calorimeter = energy_from_calorimeter + particle->energy_layer_n(0) + particle->energy_layer_n(1);
        }
        if(detected_HAD == true)
        {
          energy_detected = energy_detected + particle->energy_layer_n(2) + particle->energy_layer_n(3);
          energy_from_calorimeter = energy_from_calorimeter + particle->energy_layer_n(2) + particle->energy_layer_n(3);
        }
        if(detected_INNER_TRACKER == true)
        {
          energy_detected = energy_detected + particle->energy_layer_n(4);
          energy_from_tracker = energy_from_tracker + particle->energy_layer_n(4);
        }
        if(detected_OUTER_TRACKER == true)
        {
          energy_detected = energy_detected + particle->energy_layer_n(5);
          energy_from_tracker = energy_from_tracker + particle->energy_layer_n(5);
        }
        if(detected_STRIP_LAYER_TRACKER == true)
        {
          energy_detected = energy_detected + particle->energy_layer_n(6);
          energy_from_tracker = energy_from_tracker + particle->energy_layer_n(6);
        }
        if(detected_INNER_MUON == true)
        {
          energy_detected = energy_detected + particle->energy_layer_n(7);
          energy_from_muon_chamber = energy_from_muon_chamber + particle->energy_layer_n(7);
        }
        if(detected_OUTER_MUON)
        {
          energy_detected = energy_detected + particle->energy_layer_n(8);
          energy_from_muon_chamber = energy_from_muon_chamber + particle->energy_layer_n(8);
        }
        energies_of_all_particles_detected.push_back(energy_detected);
        energies_from_calorimeter.push_back(energy_from_calorimeter);
        energies_from_muon_chamber.push_back(energy_from_muon_chamber);
        energies_from_tracker.push_back(energy_from_tracker);
        // Finished with energy treatment
        
        if(detected_tracker == false and detected_EM == true and detected_muon_chamber == false and detected_HAD == false)
        {
          // We think this particle could be a Photon so we check against the real information and then add the particle to
          // list of correctly detected photons or incorrectly detected photon
          if(particle->get_mass() == 0 and particle->get_charge() == 0 and particle->get_spin() == 1)
          {
            identified_yet = true;
            correct_photons.push_back(particle);
            all_correct_particles.push_back(particle);
          } else
          {
            identified_yet = true;
            incorrect_photons.push_back(particle);
            all_incorrect_particles.push_back(particle);
          }
        } // Finished checking for Photon
        if(detected_tracker == true)
        {
          if(detected_EM == true and detected_HAD == false) // For Electron
          {
            // We think this particle must be an electron
            if(particle->get_mass() == 0.511)
            {
              identified_yet = true;
              correct_electrons.push_back(particle);
              all_correct_particles.push_back(particle);
            } else
            {
              identified_yet = true;
              incorrect_electrons.push_back(particle);
              all_incorrect_particles.push_back(particle);
            }
          }
          if(detected_EM == true and detected_HAD == true) // For Proton
          {
            // We think this must be a Proton
            if(particle->get_mass() == 938)
            {
              identified_yet = true;
              correct_protons.push_back(particle);
              all_correct_particles.push_back(particle);
            } else
            {
              identified_yet = true;
              incorrect_protons.push_back(particle);
              all_incorrect_particles.push_back(particle);
            }
          }
          if(detected_calorimeter == false and detected_muon_chamber == true) // For Muon
          {
            // We think this must be a muon
            if(particle->get_mass() == 105.66)
            {
              identified_yet = true;
              correct_muons.push_back(particle);
              all_correct_particles.push_back(particle);
            } else
            {
              identified_yet = true;
              incorrect_muons.push_back(particle);
              all_incorrect_particles.push_back(particle);
            }
          }

          if(detected_calorimeter == false and detected_muon_chamber == false) // Only Tracker Signature
          {
            cout<<"........................................"<<endl;
            cout<<"Only Tracker made a detection, particle can be either Muon, Electron or Proton"<<endl;
            cout<<"........................................"<<endl;
            int random{rand() % 3 + 1};
            if(random == 1)
            {
              // Electron Case
              if(particle->get_mass() == 0.511)
              {
                identified_yet = true;
                correct_electrons.push_back(particle);
                all_correct_particles.push_back(particle);
              } else
              {
                identified_yet = true;
                incorrect_electrons.push_back(particle);
                all_incorrect_particles.push_back(particle);
              }
            }
            if(random == 2)
            {
              // Muon case
              if(particle->get_mass() == 105.66)
              {
                identified_yet = true;
                correct_muons.push_back(particle);
                all_correct_particles.push_back(particle);
              } else
              {
                identified_yet = true;
                incorrect_muons.push_back(particle);
                all_incorrect_particles.push_back(particle);
              }
            }
            if(random == 3)
            {
              // Proton case
              if(particle->get_mass() == 938)
              {
                identified_yet = true;
                correct_protons.push_back(particle);
                all_correct_particles.push_back(particle);
              } else
              {
                identified_yet = true;
                incorrect_protons.push_back(particle);
                all_incorrect_particles.push_back(particle);
              }
            }
          }
        } // Finished distinguishing tracker signatures
        if(detected_tracker == false and detected_muon_chamber == false and detected_HAD == true and detected_EM == true)
        {
          // We think this must be a neutron
          if(particle->get_mass() == 939.5)
          {
            identified_yet = true;
            correct_neutrons.push_back(particle);
            all_correct_particles.push_back(particle);
          } else
          {
            identified_yet = true;
            incorrect_neutrons.push_back(particle);
            all_incorrect_particles.push_back(particle);
          }
        } // Finished looking for Neutron
        if(identified_yet == false) // IF PARTICLE STILL NOT IDENTIFIED THEN IT GOES TO GHOST PARTICLE VECTOR
        {
          if(particle->get_mass() == 0 and particle->get_spin() == 0.5)
          {
            neutrinos_through_detector.push_back(particle);
            hidden_energy.push_back(particle->get_true_energy());
          } else
          {
            ghost_particles.push_back(particle);
            hidden_energy.push_back(particle->get_true_energy());
          }
        }
      } // Finished iterating particles in the event

      // Sumations
      double sum_of_energies{0};
      double sum_of_hidden_energies{0};
      for(double num : energies_of_all_particles_detected)
      {
        sum_of_energies = sum_of_energies + num;
      }
      for(double num : hidden_energy)
      {
        sum_of_hidden_energies = sum_of_hidden_energies + num;
      }
      num_of_correct_electrons = num_of_correct_electrons + correct_electrons.size();
      num_of_correct_muons = num_of_correct_muons + correct_muons.size();
      num_of_correct_protons = num_of_correct_protons + correct_protons.size();
      num_of_correct_photons = num_of_correct_photons + correct_photons.size();
      num_of_correct_neutrons = num_of_correct_neutrons + correct_photons.size();
      num_of_incorrect_electrons = num_of_incorrect_electrons + incorrect_electrons.size();
      num_of_incorrect_muons = num_of_incorrect_muons + incorrect_muons.size();
      num_of_incorrect_protons = num_of_incorrect_protons + incorrect_protons.size();
      num_of_incorrect_photons = num_of_incorrect_photons + incorrect_photons.size();
      num_of_incorrect_neutrons = num_of_incorrect_neutrons + incorrect_photons.size();
      num_of_neutrinos = num_of_neutrinos + neutrinos_through_detector.size();
      num_of_ghosts = num_of_ghosts + ghost_particles.size();

      // PRINTING EVERYTHING AFTER THE EVENT
      cout<<"//////////////////////////"<<endl;
      cout<<"Event Information:"<<events_analysed<<endl;
      cout<<".................."<<endl;
      cout<<"Particles Detected: "<<total_particles_detected<<endl;
      cout<<"-------------------"<<endl;
      int count{0};
      for(std::shared_ptr<Particle> elem : all_particles_detected)
      {
        cout<<"Particle "<<elem->get_identifier()<<endl;
        elem->print_name();
        cout<<"Energy from Tracker: "<<energies_from_tracker[count]<<endl;
        cout<<"Energy from Calorimeter: "<<energies_from_calorimeter[count]<<endl;
        cout<<"Energy from Muon Chamber: "<<energies_from_muon_chamber[count]<<endl;
        count++;
        cout<<"---------------"<<endl;
      }

      cout<<"Overall for Event: "<<events_analysed<<endl;
      cout<<"__________________"<<endl;
      cout<<"The total energy detected was: "<<sum_of_energies<<endl;
      cout<<"The total energy missed was: "<<sum_of_hidden_energies<<endl;
      cout<<"Percentage of energy detected: "<<((float)sum_of_energies/(float)(sum_of_hidden_energies + sum_of_energies))*100<<"%"<<endl;
      cout<<"Efficiency of run (detected particles/event particles): "<<((float)total_particles_detected / (float)number_of_particles_in_event) * 100<<"%"<<endl;
      cout<<"True efficiency of run: (correctly detected / event particles): " << ((float)(total_particles_detected - all_incorrect_particles.size()) / (float)number_of_particles_in_event) * 100<<"%"<<endl;
      cout<<"Correctly Detected / Detected: "<<((float)(total_particles_detected - all_incorrect_particles.size()) / (float)total_particles_detected) * 100<<"%"<<endl;
      cout<<"Incorrect Protons: "<<incorrect_protons.size()<<endl;
      cout<<"Incorrect Photons: "<<incorrect_photons.size()<<endl;
      cout<<"Incorrect Electrons: "<<incorrect_electrons.size()<<endl;
      cout<<"Incorrect Muons: "<<incorrect_muons.size()<<endl;
      cout<<"Incorrect Neutrons: "<<incorrect_neutrons.size()<<endl;
      cout<<"//////////////////////////"<<endl; cout<<endl;
      
      // Recording for all the events that passed through
      particles_per_event.push_back(all_particles_detected); 
      correct_particles_per_event.push_back(all_correct_particles); 
      incorrect_particles_per_event.push_back(all_incorrect_particles); 
      neutrinos_per_event.push_back(neutrinos_through_detector);
      ghosts_per_event.push_back(ghost_particles);
      energies_of_all_per_event.push_back(energies_of_all_particles_detected);
      hidden_energy_per_event.push_back(hidden_energy);
      energies_from_tracker_per_event.push_back(energies_from_tracker);
      energies_from_calorimeter_per_event.push_back(energies_from_calorimeter);
      energies_from_muon_chamber_per_event.push_back(energies_from_muon_chamber);
      events_analysed++; // STATIC
      //
      this->clear();
      return true;
    } 

    // This function will set everything to the default value.
    void clear()
    {
      has_detected = false;
      for (int i{0}; i < trackers_in_detector.size(); i++)
      {
        trackers_in_detector[i]->clear_tracker();
      }
      for (int i{0}; i < calorimeters_in_detector.size(); i++)
      {
        calorimeters_in_detector[i]->clear_calorimeter();
      }
      for (int i{0}; i < muon_chambers_in_detector.size(); i++)
      {
        muon_chambers_in_detector[i]->clear_muon_chamber();
      }

      hidden_energy.clear();
      ghost_particles.clear();
      neutrinos_through_detector.clear();
      correct_muons.clear();
      correct_protons.clear();
      correct_neutrons.clear();
      correct_photons.clear();
      correct_electrons.clear();
      incorrect_muons.clear();
      incorrect_protons.clear();
      incorrect_neutrons.clear();
      incorrect_photons.clear();
      incorrect_electrons.clear();
      energies_of_all_particles_detected.clear();
      all_correct_particles.clear();
      all_incorrect_particles.clear();
      total_particles_detected = 0;
      all_particles_detected.clear();
      number_of_particles_in_event = 0;
    }
};

// Useful Functions
void distribute_energy(vector<std::shared_ptr<Particle>> event, double energy) // This function helps distributing the initial energy of the event randomly
{
  // Make exception handling for the energy...

  std::vector<int> divisions(event.size()); // vector in which we store fractions of the energy
  for (int i{0}; i < event.size() - 1; ++i)
  {
    divisions[i] = 1 + rand() % ((int)energy - event.size() + 1); 
    energy = energy - divisions[i]; 
  }
  divisions.back() = energy;

  int count{0};
  for(std::shared_ptr<Particle> particle : event)
  {
    particle->set_true_energy(divisions[count]);
    count = count + 1;
  }

}

int main()
{
  srand(time(NULL)); // For changing seed overall (we also have seed changing when running to completely avoid repetition)
  seed = rand();

  Detector D1;
  vector<std::shared_ptr<Particle>> event;
  vector<std::shared_ptr<Particle>> event2;
  vector<std::shared_ptr<Particle>> event3;
  vector<std::shared_ptr<Particle>> event4;
/*
  // EVENT NUMBER 1: Two Photons
  // ..............
  std::shared_ptr<Photon> ph1_ptr = std::make_shared<Photon>();
  std::shared_ptr<Photon> ph2_ptr = std::make_shared<Photon>();
  event.push_back(ph1_ptr);
  event.push_back(ph2_ptr);
  D1.detect(event);
*/
/*
  // EVENT NUMBER 2: Proton - Neutron - Photon
  // ..............
  std::shared_ptr<Neutron> n1_ptr = std::make_shared<Neutron>();
  std::shared_ptr<Photon> ph3_ptr = std::make_shared<Photon>();
  std::shared_ptr<Proton> pr1_ptr = std::make_shared<Proton>();
  event2.push_back(n1_ptr);
  event2.push_back(ph3_ptr);
  event2.push_back(pr1_ptr);
  D1.detect(event2);
*/
  // EVENT NUMBER 3: Proton - Neutron - Photon // AQUI HAY ERRORES// AQUI HAY ERRORES// AQUI HAY ERRORES// AQUI HAY ERRORES
  // ..............
  std::shared_ptr<Tau> t1_ptr = std::make_shared<Tau>();
  std::shared_ptr<Tau> t2_ptr = std::make_shared<Tau>();
  event3.push_back(t1_ptr);
  event3.push_back(t2_ptr);
  D1.detect(event3);
  if(t1_ptr->get_decay_type() == "leptonic")
  {
    for(int i{0}; i < 3; i++)
    {
      t1_ptr->get_decay_products()[i]->print_name();
    }
  }
  if(t2_ptr->get_decay_type() == "leptonic")
  {
    for(int i{0}; i < 3; i++)
    {
      t2_ptr->get_decay_products()[i]->print_name();
    }
  }
/*
  // EVENT NUMBER 4: Proton - Neutron - Photon
  // ..............
  std::shared_ptr<Proton> pr2_ptr = std::make_shared<Proton>();
  std::shared_ptr<Neutrino> neutri1_ptr = std::make_shared<Neutrino>();
  std::shared_ptr<Neutrino> neutri2_ptr = std::make_shared<Neutrino>();
  event4.push_back(n1_ptr);
  event4.push_back(ph3_ptr);
  event4.push_back(pr1_ptr);
  D1.detect(event4);
*/
  return 0;
}