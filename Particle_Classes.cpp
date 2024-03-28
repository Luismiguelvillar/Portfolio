/*
Assignment 4:
Enhancing particle class to describe some
leptons (electrons, muons, antielectrons and antimuons) of
the standard model of physics.

by Luis Miguel Villar Padruno
3/3/2024
*/

#include <iostream>
#include <string>
#include <vector>
#include<math.h>

using std::cout;
using std::endl;
using std::cin;
using std::string;
using std::vector;

//Please Ignore this (preparation for next assignment)

class FourMomentum 
{
  private:
    double Energy_Over_c{0};
    double x1_coordinate{0};
    double x2_coordinate{0};
    double x3_coordinate{0};
  public:
    FourMomentum() = default;
    FourMomentum(double Energy_Over_c_initial, double x1_initial, double x2_initial, double x3_initial)
    {
      if(Energy_Over_c_initial < 0)
      {
        cout << "Invalid value for the Energy/c: ";
        cin>>Energy_Over_c_initial;
        while (cin.fail() or Energy_Over_c_initial < 0)
        {
          cout<<"Invalid Energy/c, try again:  ";
          cin>>Energy_Over_c_initial;
        }
        Energy_Over_c = Energy_Over_c_initial;
      } else
      {
        Energy_Over_c = Energy_Over_c_initial;
      }
      x1_coordinate = x1_initial;
      x2_coordinate = x2_initial;
      x3_coordinate = x3_initial;
    }
    ~FourMomentum(){cout<<"Four Momentum vector destroyed"<<endl;}
    void print_FourMomentum();
};

void FourMomentum::print_FourMomentum()
{
    cout<<"Four Momentum = ("<<Energy_Over_c<<", "<<x1_coordinate<<", "<<x2_coordinate<<", "<<x3_coordinate<<")"<<endl;
}

class Particle 
{
  private:
    const double light_speed{2.99792458e8};
    double mass {0.0};
    double charge {0.0};
    string name {"ghost"};
    vector<double>* four_momentum{nullptr};
    string custom_setter{"default"}; // default  or custom
  public:
    Particle()
    {
      cout<<"Default Constructor Called"<<endl;
      four_momentum = new vector<double> (4);
      const double light_speed{2.99792458e8};
      double mass {0.0};
      double charge {0.0};
      string name {"ghost"};
      vector<double>* four_momentum{nullptr};
      string custom_setter{"default"}; // default  or custom
    }
    Particle(string type_initial, double mass_initial, double charge_initial, double energy_over_c, double x1, double x2, double x3, string custom_setter_initial)
    {
      cout<<"Parametrized Constructor Called."<<endl;
      if(custom_setter_initial != "default" and custom_setter_initial != "custom")
      {
        cout<<"Custom setter must be true or false: ";
        cin>>custom_setter_initial;
        while(cin.fail() or (custom_setter_initial != "default" and custom_setter_initial != "custom"))
        {
          cout<<"Invalid input, try again: ";
          cin.clear();
          cin.ignore();
          cin>>custom_setter_initial;
        }
        custom_setter = custom_setter_initial;
      } else
      {
        custom_setter = custom_setter_initial;
      }
      if(type_initial != "muon" and type_initial != "electron" and type_initial != "antielectron" and type_initial != "antimuon")
      {
        cout<<"invalid name for the particle, try again: ";
        cin>>type_initial;
        while(cin.fail() or (type_initial != "muon" and type_initial != "electron" and type_initial != "antielectron" and type_initial != "antimuon"))
        {
          cout << "invalid name, try again: ";
          cin.clear();
          cin.ignore();
          cin>>type_initial;
        }
        name = type_initial;
        
      } else
      {
        name = type_initial;
      }
      if(mass_initial < 0)
      {
        cout<<"invalid mass for particle, try again: ";
        cin>>mass_initial;
        while(cin.fail() or (mass_initial < 0))
        {
          cout << "invalid mass, try again: ";
          cin.clear();
          cin.ignore();
          cin>>mass_initial;
        }
        if(mass_initial != 0.511 and mass_initial != 105.66)
        {
          string selection{""};
          cout<<"The particle you input is a/an "<<type_initial<<" but the mass does not match the standard model.";
          cout<<" Do you want to keep your selection or use the correct value for the mass? (correct/selection): ";
          cin>>selection;
          while(cin.fail() or (selection != "correct" and selection != "selection"))
          {
            cout<<"Invalid input";
            cin.clear();
            cin.ignore();
            cin>>selection;
          }
          if(selection == "correct")
          {
            if(type_initial == "electron" or type_initial == "antielectron")
            {
              mass = 0.511;
            } else
            {
              mass = 105.66;
            }
          } else
          {
            mass = mass_initial;
          }
          
        } else
        {
          mass = mass_initial;
        }
        
      } else
      {
        mass = mass_initial;
      }

      if(charge_initial != 1 and charge_initial != -1)
      {
        cout<<"Invalid charge, try again: ";
        cin>>charge_initial;
        while(cin.fail() or (charge_initial != 1 and charge_initial != -1))
        {
            cout<<"Invalid charge, try again: ";
            cin.clear();
            cin.ignore();
            cin>>charge_initial;
        }
        if(((type_initial == "electron" or type_initial == "muon") and (charge_initial == -1)) or ((type_initial == "antielectron" or type_initial == "antimuon") and (charge_initial == 1)))
        {
          charge = charge_initial;
        } else
        {
          string selection{""};
        cout<<"The particle you input is a/an "<<type_initial<<" but the charge does not match the standard model.";
        cout<<" Do you want to keep your selection or use the correct value for the charge? (correct/selection): ";
        cin>>selection;
        while(cin.fail() or (selection != "correct" and selection != "selection"))
        {
        cout<<"Invalid input (correct/selection): ";
        cin.clear();
        cin.ignore();
        cin>>selection;
        }
        if(selection == "correct")
        {
          if(type_initial == "electron" or type_initial == "muon")
          {
            charge = -1;
          } else
          {
            charge = 1;
          }
        } else
        {
          charge = charge_initial;
        }
        }
      } else
      {
        charge = charge_initial;
      }
      // Four Momentum treatment 
      four_momentum = new vector<double> (4);
      four_momentum->at(1) = x1;
      four_momentum->at(2) = x2;
      four_momentum->at(3) = x3;
      if(energy_over_c > 0)
      {
        four_momentum->at(0) = energy_over_c;
      } else
      {
        cout<<"No negative energy allowed: ";
        cin>>energy_over_c;
        while(cin.fail() or energy_over_c < 0)
        {
            cout<<"Invalid value, try again: ";
            cin.clear();
            cin.ignore();
            cin>>energy_over_c;
        }
        four_momentum->at(0) = energy_over_c;
      }
    }
    ~Particle(){ cout<<"Destructor called. "; cout<<"Particle destroyed and memory freed."<<endl; delete four_momentum;}
    //Print name class
    void print_name();
    //Set functions
    void set_name(const string name_set)
    {
      if(name_set == "electron" or name_set == "muon" or name_set == "antimuon" or name_set == "antielectron")
      {
        name = name_set;
        if((name_set == "antimuon" or name_set == "antielectron") and custom_setter == "default")
        {
          charge = 1;
        } else 
        {
          charge = -1;
        }
        if((name_set == "muon" or name_set == "antimuon") and custom_setter == "default")
        {
          mass = 105.66;
        } else
        {
          mass = 0.511;
        }
      }
      else
      {
        cout << "The lepton number was not recognized" << endl;
        name = name_set;   
      }
      return;
    }
    void set_charge(const double charge_set) 
    {
      if(charge_set == 1 or charge_set == -1)
      {
        charge = charge_set;
        if(charge_set == 1 and (name != "antielectron" and name != "antimuon") and custom_setter == "default") 
        {
          string anti {"anti"};
          name = anti.append(name);
        } 
        if(charge_set == -1 and (name == "antielectron" or name == "antimuon") and custom_setter == "default")
        {
          name.erase(0,4);
        }
      } else 
      {
        if(custom_setter == "default")
        {
          cout << "The charge of the lepton should be 1 or -1" << endl;
          if(name == "electron" or name == "muon")
          {
            charge = -1;
          } else if(name == "antimuon" or name == "antielectron")
          {
            charge = 1;
          } else
          {
            charge = 0;
          }
        } else
        {
          charge = charge_set;
        }
      }
      return;
    }
    void set_mass(const double mass_set) 
    {
      if(mass >= 0)
      {
        mass = mass_set;
        if(custom_setter == "default")
        {
          if((std::abs(mass_set- 0.511) < std::abs(mass_set - 105.66)))
          {
            if(charge == -1)
            {
              name = "electron";
            } else
            {
              name = "antielectron";
            }
          } else
          {
            if(charge == -1)
            {
              name = "muon";
            } else
            {
              name = "antimuon";
            }
          }
        } else
        {
          mass = mass_set;
        }
      } else 
      {
        cout << "Mass of Particle should be positive" << endl;
        mass = 0;
      }
      return;
    }
    void set_custom_setter(const string custom_setter_set)
    {
      if(custom_setter_set == "default" or custom_setter_set == "custom")
      {
        custom_setter = custom_setter_set;
      } else
      {
        custom_setter = "default";
      }
    }
    void set_energy_over_c(const double energy_set)
    {
      if(custom_setter != "default")
      {
        four_momentum->at(0) = energy_set;
        return;
      } else
      {
        if(energy_set < 0)
        {
            cout<<"Invalid value for energy set, the particle with adquire the absolute value of the given value."<<endl;
            four_momentum->at(0) = fabs(energy_set);
            return;
        } else
        {
         four_momentum->at(0) = energy_set;
         return;
        }
      }
      return;
    }
    void set_momentum_px(const double px_set)
    {
      four_momentum->at(1) = px_set;
      return;
    }
    void set_momentum_py(const double py_set)
    {
      four_momentum->at(2) = py_set;
      return;
    }
    void set_momentum_pz(const double pz_set)
    {
      four_momentum->at(3) = pz_set;
      return;
    }
    //Get functions
    string get_name() const {return name;}
    double get_charge() const {return charge;}
    double get_mass() const {return mass;}
    double get_energy_over_c() const {return four_momentum->at(0);}
    double get_momentum_px() const {return four_momentum->at(1);}
    double get_momentum_py() const {return four_momentum->at(2);}
    double get_momentum_pz() const {return four_momentum->at(3);}

    // Overloading operators
    Particle operator+(Particle const & object)
    {
      Particle temp{name, mass, charge, four_momentum->at(0), four_momentum->at(1), four_momentum->at(2), four_momentum->at(3), custom_setter};
      temp.set_energy_over_c(four_momentum->at(0) +  object.four_momentum->at(0));
      temp.set_momentum_px(four_momentum->at(1) + object.four_momentum->at(1));
      temp.set_momentum_py(four_momentum->at(2) + object.four_momentum->at(2));
      temp.set_momentum_pz(four_momentum->at(3) + object.four_momentum->at(3));
      return temp; // returns a new particle of the same type as the dirst one but with the sum of the four momenta
    }

    // Dot product and JUST four-momentum summation
    double dotProduct(Particle const & object) 
    {
      return four_momentum->at(0) * object.four_momentum->at(0) - four_momentum->at(1) * object.four_momentum->at(1) - four_momentum->at(2) * object.four_momentum->at(2) - four_momentum->at(3) * object.four_momentum->at(3); 
    }
    vector<double> sum_four_vector(Particle const & object)
    {
      vector<double> temporal_vector (4);
      temporal_vector[0] = four_momentum->at(0) + object.four_momentum->at(0);
      temporal_vector[1] = four_momentum->at(1) + object.four_momentum->at(1);
      temporal_vector[2] = four_momentum->at(2) + object.four_momentum->at(2);
      temporal_vector[3] = four_momentum->at(3) + object.four_momentum->at(3);
      return temporal_vector; 
    }

    //Visualize the four vector
    void print_four_vector()
    {
      cout<<"["<<four_momentum->at(0)<<", "<<four_momentum->at(1)<<", "<<four_momentum->at(2)<<", "<<four_momentum->at(3)<<"]"<<endl;
      return;
    }

    // Costumize Special Functions
    Particle & operator=(Particle const & object)
    {
      if(&object == this) return *this;
      cout<<"Calling Assignment Operator."<<endl;
      delete four_momentum;
      four_momentum = new vector<double> (4);
      name = object.name;
      mass =  object.mass;
      charge = object.charge;
      four_momentum->at(0) = object.four_momentum->at(0);
      four_momentum->at(1) = object.four_momentum->at(1);
      four_momentum->at(2) = object.four_momentum->at(2);
      four_momentum->at(3) = object.four_momentum->at(3);
      custom_setter = object.custom_setter;
      return *this;
    } 
    Particle(Particle const & object)
    {
      cout<<"Calling Copy Constructor."<<endl;
      four_momentum = new vector<double> (4);
      name = object.name;
      mass =  object.mass;
      charge = object.charge;
      four_momentum->at(0) = object.four_momentum->at(0);
      four_momentum->at(1) = object.four_momentum->at(1);
      four_momentum->at(2) = object.four_momentum->at(2);
      four_momentum->at(3) = object.four_momentum->at(3);
      custom_setter = object.custom_setter;
    }
    Particle & operator=(Particle && particle_to_migrate)
    {
      cout<<"Calling Move assignment operator."<<endl;
      std::swap(name, particle_to_migrate.name);
      std::swap(mass, particle_to_migrate.mass);
      std::swap(charge, particle_to_migrate.charge);
      std::swap(four_momentum->at(0), particle_to_migrate.four_momentum->at(0));
      std::swap(four_momentum->at(1), particle_to_migrate.four_momentum->at(1));
      std::swap(four_momentum->at(2), particle_to_migrate.four_momentum->at(2));
      std::swap(four_momentum->at(3), particle_to_migrate.four_momentum->at(3));
      std::swap(custom_setter, particle_to_migrate.custom_setter);
      return *this;

    }
    Particle(Particle && object)
    {
      cout<<"Calling Move Constructor."<<endl;
      four_momentum = new vector<double> (4);
      name = object.name;
      mass =  object.mass;
      charge = object.charge;
      four_momentum->at(0) = object.four_momentum->at(0);
      four_momentum->at(1) = object.four_momentum->at(1);
      four_momentum->at(2) = object.four_momentum->at(2);
      four_momentum->at(3) = object.four_momentum->at(3);
      custom_setter = object.custom_setter;
      //Now we clear the memory of the moved particle
      object.name = "ghost";
      object.mass = 0;
      object.charge = 0;
      object.four_momentum->at(0) = 0;
      object.four_momentum->at(1) = 0;
      object.four_momentum->at(2) = 0;
      object.four_momentum->at(3) = 0;
      object.custom_setter = "default";
    }
}; // end of particle class

void Particle::print_name() 
{
  cout<<"Particle: [type, mass, charge, E/c, px, py, pz, custom] = "<<"["<<name<<", "<<mass<<", "<<charge<<", "<<four_momentum->at(0)<<", "<<four_momentum->at(1)<<", "<<four_momentum->at(2)<<", "<<four_momentum->at(3)<<", "<<custom_setter<< "]"<<endl;
  return;
}

class detector
{
  private:
    string detector_type{"tracker"};
    bool is_on{true};
    int number_of_particles_detected{0};
    int number_of_muons{0};
    int number_of_electrons{0};
    int number_of_antielectrons{0};
    int number_of_antimuons{0};

  public:
  detector() = default;
  detector(string initial_detector_type, bool initial_is_on, int initial_particles_detected): 
  detector_type{initial_detector_type}, is_on{initial_is_on} {}
  ~detector() {cout<<"Detector destroyed"<<endl;}
  // Print function
  void print_detector()
  {
    cout<<"Detector: [type, status, particles detected, electrons, antielectrons, muons, antimuons] = "<<"["<<detector_type<<", "<<is_on<<", "<<number_of_particles_detected<<", "<<number_of_electrons<<", "<<number_of_antielectrons<<", "<<number_of_muons<<", "<<number_of_antimuons<<"]"<<endl;
  }
  // Set Functions
  void set_detector_type(string const detector_type_set)
  {
    if(detector_type_set == "tracker" or detector_type_set == "calorimeter" or detector_type_set == "muon chamber")
    {
      detector_type = detector_type_set;
    } else
    {
      cout << "detector type not recognized, default is tracker" << endl;
      detector_type = "tracker";
    }
    return;
  }
  void set_is_on(bool const is_on_set)
  {
    if(is_on_set == true or is_on_set == false)
    {
      is_on = is_on_set;
    } else 
    {
      cout << "Error setting status of detector, turning off" << endl;
      is_on = false;
    }
    return;
  }
  void clean()
  {
    number_of_antielectrons = 0;
    number_of_antimuons = 0;
    number_of_electrons = 0;
    number_of_muons = 0;
    number_of_particles_detected = 0;
    return;
  }
  // Detection Process
  bool detection(const Particle& p) 
  {
    if(is_on == true)
    {
      if(detector_type == "tracker") // All
      {
        if(p.get_name() != "ghost")
        {
          if(p.get_charge() == -1) 
          {
            if(p.get_name() == "electron")
            {
              number_of_electrons =  number_of_electrons + 1;
              number_of_particles_detected = number_of_particles_detected + 1;
              cout<<"electron detected. (Via tracker)"<<endl; 
              return true;
            } else
            {
              number_of_muons = number_of_muons + 1;
              number_of_particles_detected = number_of_particles_detected + 1;
              cout<<"muon detected. (Via tracker)"<< endl;
              return true;
            }
          } else
          {
            if(p.get_name() == "antielectron")
            {
              number_of_antielectrons =  number_of_antielectrons + 1;
              number_of_particles_detected = number_of_particles_detected + 1;
              cout<<"antielectron detected. (Via tracker)"<<endl;
              return true;
            } else
            {
              number_of_antimuons = number_of_antimuons + 1;
              number_of_particles_detected = number_of_particles_detected + 1;
              cout<<"antimuon detected. (Via tracker)"<<endl;
              return true;
            }
          }
        } else 
        {
          return false;
        }
      }
      if(detector_type == "calorimeter") // electrons & antielectrons
      {
        if(p.get_name() != "ghost")
        {
          if(p.get_name() == "electron" and p.get_charge() == -1)
          {
            number_of_electrons =  number_of_electrons + 1;
            number_of_particles_detected = number_of_particles_detected + 1;
            cout<<"electron detected. (Via calorimeter)"<<endl;
            return true;
          } else 
          {
            if(p.get_name() == "antielectron" and p.get_charge() == 1)
            {
              number_of_antielectrons =  number_of_antielectrons + 1;
              number_of_particles_detected = number_of_particles_detected + 1;
              cout<<"antielectron detected. (Via calorimeter)"<<endl;
              return true;
            } else
            {
              return false;
            }
          }
        } else
        {
          return false;
        }
      }
      if(detector_type == "muon chamber") // muons & antimuons
      {
        if(p.get_name() != "ghost")
        {
          if(p.get_name() == "muon" and p.get_charge() == -1) 
          {
            number_of_muons = number_of_muons + 1;
            number_of_particles_detected = number_of_particles_detected + 1;
            cout<<"muon detected. (Via muon chamber)"<<endl;
            return true;
          } else 
          {
            if(p.get_name() == "antimuon" and p.get_charge() == 1)
            {
              number_of_antimuons = number_of_antimuons + 1;
              number_of_particles_detected = number_of_particles_detected + 1;
              cout<<"antimuon detected. (Via muon chamber)"<<endl;
              return true;
            } else
            {
              return false;
            }
          }
        } else
        {
          return false;
        }
      } else{
        return false;
      }
    } else 
    {
      return false;
    }
  }
  // Get functions
  int get_num_of_particles_detected() const {return number_of_particles_detected;}
  int get_num_of_muons() const {return number_of_muons;}
  int get_num_of_electrons() const {return number_of_electrons;}
  int get_num_of_antielectrons() const {return number_of_antielectrons;}
  int get_num_of_antimuons() const {return number_of_antimuons;}
  bool get_status() const {return is_on;}
};

int main ()
{
  Particle p;
  vector<Particle> test_vector;

  // Generating vector of particles (e, e, m, m, m, m, ae, am)
  for(int i{0}; i < 8; i++)
  {
    test_vector.push_back(p);
    if(i == 0 or i == 1)
    {
      test_vector[i].set_name("electron");
      test_vector[i].set_energy_over_c(10+5*i);
      test_vector[i].set_momentum_px(1+2*i);
      test_vector[i].set_momentum_py(3-i);
      test_vector[i].set_momentum_pz(5+3*i);
    } else if(i <= 5 and i > 1)
    {
      test_vector[i].set_name("muon");
      test_vector[i].set_energy_over_c(10+5*i);
      test_vector[i].set_momentum_px(1+2*i);
      test_vector[i].set_momentum_py(3-i);
      test_vector[i].set_momentum_pz(5+3*i);
    } else if(i == 6)
    {
      test_vector[i].set_name("antielectron");
      test_vector[i].set_energy_over_c(10+5*i);
      test_vector[i].set_momentum_px(1+2*i);
      test_vector[i].set_momentum_py(3-i);
      test_vector[i].set_momentum_pz(5+3*i);
    } else if(i == 7)
    {
      test_vector[i].set_name("antimuon");
      test_vector[i].set_energy_over_c(10+5*i);
      test_vector[i].set_momentum_px(1+2*i);
      test_vector[i].set_momentum_py(3-i);
      test_vector[i].set_momentum_pz(5+3*i);
    }
  }


  //Visualize particles
  cout<<endl;
  cout<<"--- Particles we have created in the vector: ---"<<endl;
  for(int j{0}; j < test_vector.size(); j++)
  {
    test_vector[j].print_name();
  }

  // Now we test...
  // --------------
  //Summing two four-momenta of electrons
  cout<<endl;
  cout<<"The four momentum of each electron is:"<<endl;
  cout<<"Electron 1: ";
  test_vector[0].print_four_vector();
  cout<<"Electron 2: ";
  test_vector[1].print_four_vector();
  cout<<"The sum of the vectors is: "<<endl;
  (test_vector[0] + test_vector[1]).print_four_vector();

  //Dot Product of first two muons [Note that the dot product is defined as p*p = (E/c)^2 - px^2 - py^2 - pz^2]
  cout<<"The four momentum of each muon is:"<<endl;
  cout<<"Muon 1: ";
  test_vector[2].print_four_vector();
  cout<<"Muon 2: ";
  test_vector[3].print_four_vector();
  cout<<"The dot product is: " << endl;
  cout<<test_vector[2].dotProduct(test_vector[3])<<endl;
  cout<<endl;

  //Assignment operator of an electron to a new electron
  Particle p_new;
  cout<<"Before assigning p_new is: "; 
  p_new.print_name();
  cout<<"and the copied electron is: ";
  test_vector[0].print_name();
  p_new = test_vector[0];
  cout<<"After assignment operator: ";
  p_new.print_name();
  cout<<"and the copied electron is: ";
  test_vector[0].print_name(); 
  cout<<"-TEST OF ASSIGNING-"<<endl;
  cout<< "If we change values in the original electron, no changes will be seen in the p_new particle."<<endl;
  test_vector[0].set_momentum_py(100);
  cout<<"The py of the copied electron was set to 100: "<< endl;
  cout<<"Original electron: ";
  test_vector[0].print_name();
  cout<<"p_new electron: ";
  p_new.print_name();
  cout<<endl;

  // Copy constructor of the first muon to a new muon
  cout<<"Now we will create a new muon using the Copy Constructor: "<<endl;
  cout<<"We eill copy the following muon: ";
  test_vector[2].print_name();
  Particle new_muon{test_vector[2]};
  cout<<"The new Particle is: ";
  new_muon.print_name();
  cout<<"-TEST OF COPY OPERATOR-"<<endl;  
  cout<<"If we change the values of the original muon, no changes will be seen in new_muon."<<endl;
  test_vector[2].set_energy_over_c(1000);
  cout<<"The E/c of the original muon was set to 1000. "<< endl;
  cout<<"Original muon: ";
  test_vector[2].print_name();
  cout<<"New muon: ";
  new_muon.print_name();
  cout<<endl;

  // Move antielectron into another antielectron using Move Constructor
  cout<<"Now we use the Move Constructor to Move the antielectron to another antielectron"<<endl;
  cout<<"This operation will leave the original antielectron empty."<<endl;
  cout<<"Antielectron before Move Constructor: ";
  test_vector[6].print_name();
  Particle new_antielectron {std::move(test_vector[6])};
  cout<<"New antielectron after Move Constructor: ";
  new_antielectron.print_name();
  cout<<"Original Antielectron after the Move Constructor: ";
  test_vector[6].print_name();
  cout<<endl;

  // Assing antimuon to another antimuon using Move Assignment
  cout<<"Now we assign the antimuon to another using Move Assignemnt: "<<endl;
  cout<<"The antimuon before the Move Assignment: ";
  test_vector[7].print_name();
  Particle new_antimuon;
  cout<<"The new antimuon before Move Assignment: ";
  new_antimuon.print_name();
  new_antimuon = std::move(test_vector[7]);
  cout<<"Original antimuon after Move Assignment: ";
  test_vector[7].print_name();
  cout<<"New antimuon after Move Assignment: ";
  new_antimuon.print_name();



  cout<<endl;
  //Clear Memory
  test_vector.clear();
  return 0;
}
