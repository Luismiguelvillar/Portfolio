/*
Assignment 2 (PHYS30762)
------------------------
This program reads a .dat document containing information about different university
courses. The data read will then be used to calculate the mean of the grades, the standard
deviation and the standard error in the mean for a selected year (selected by the user). 
The program will output the full list of courses in alphabetical or code order, as the
user prefers.

Created by Luis Miguel Villar Padruno
15/02/2024
4:44 pm
*/

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<string>
#include<vector>
#include<typeinfo>
#include<bits/stdc++.h>

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std:: sort;

// Functions to compute mean, standard deviation or for other tasks.
float mean_calculation(vector <float> course_grades)
{
  float summation{0};

  for(int h{0}; h < course_grades.size(); h++)
  {
    summation = summation + course_grades[h];
  }
  return summation/course_grades.size();
}

float standard_deviation(vector <float> course_grades)
{
  float summation{0};
  float mean_of_grades{mean_calculation(course_grades)};
  for(int l{0}; l < course_grades.size(); l++)
  {
    summation = summation + pow(course_grades[l] - mean_of_grades, 2.0);
    
  }

  return sqrt(summation /(course_grades.size() - 1));
}

float standard_deviation_error(float std_deviation, float number_of_courses)
{
  return std_deviation / number_of_courses;
}


// Main function

int main()
{
  // Define variables
  string data_file;
  vector<string> all_courses_and_data;
  vector<float> grades;
  vector<int> course_code;
  vector<string> course_names;
  vector<char> auxiliar;
  vector<int> auxiliar_year;
  vector<string> sub_vector_all_courses_and_data;
  vector<int> sub_courses_codes;
  vector<float> sub_grades;
  vector<string> sub_courses_names;
  vector<int> indexes_code_sort;
  vector<string> ordered_all;
  vector<int> indexes_name_sort;

  vector<string> all_courses_and_data_copy;
  vector<float> grades_copy;
  vector<int> course_code_copy;
  vector<string> course_names_copy;

  int* point_cont;
  point_cont = new int;
  *point_cont = 0;
  float* course_number_pointer;
  course_number_pointer = new float;
  *course_number_pointer = 0;
  char* pointer_yes_no;
  pointer_yes_no = new char;
  *pointer_yes_no = 'y';
  float* std_deviation_pointer;
  std_deviation_pointer = new float;
  *std_deviation_pointer = 0;
  float* error_std_dev_pointer;
  error_std_dev_pointer = new float;
  *error_std_dev_pointer = 0;
  float* mean_of_grades_pointer;
  mean_of_grades_pointer = new float;
  *mean_of_grades_pointer = 0;
  int* year_selection_pointer;
  year_selection_pointer = new int;
  *year_selection_pointer = 0;
  char* alphabet_or_code_pointer;
  alphabet_or_code_pointer = new char;
  *alphabet_or_code_pointer = 'c';


  // Ask user to enter filename
  cout<<"Enter data filename: ";
  cin>>data_file;

  // Open file (you must check if successful)
  std::fstream courses(data_file);

  if (!courses.good())
  {
    std::cerr<<"Error: file could not be opened..."<<endl;
    return(1);
  }

  string line;
  while(std::getline(courses, line)) 
  {
    all_courses_and_data.push_back(line);
    *point_cont = 0;

    for (int i = 0; i < line.size(); i++)
    {
      if(line[i] == ' ') 
      { 
        if(*point_cont == 0) 
        {
          string str(auxiliar.begin(), auxiliar.end());
          grades.push_back(std::stof(str));
          auxiliar.clear();


        } else if(*point_cont == 1)
        {
          string str(auxiliar.begin(), auxiliar.end());
          course_code.push_back(std::stoi(str));
          auxiliar.clear();
          *point_cont = *point_cont+1;
        }
        *point_cont = *point_cont+1;
      }
      
      // We use the auxiliar list
      auxiliar.push_back(line[i]);
    }
    string str(auxiliar.begin(), auxiliar.end());
    str.erase(str.begin(),str.begin() + 3);
    course_names.push_back(str);
    auxiliar.clear();

  }

  // Close file 
  courses.close();

  while(*pointer_yes_no == 'y')
  {
    //Calculations and assignment
    *mean_of_grades_pointer = mean_calculation(grades);
    *std_deviation_pointer = standard_deviation(grades);
    *error_std_dev_pointer = standard_deviation_error(*std_deviation_pointer, grades.size());
    
    //Sorting order question to user
    cout<<"Are you interested in courses from a specific year? (0 / 1 / 2 / 3 / 4) (Choose 0 to see all courses): ";
    cin>> *year_selection_pointer;
    while (cin.fail() or (*year_selection_pointer != 0 and *year_selection_pointer != 1
    and *year_selection_pointer != 2 and *year_selection_pointer != 3 and *year_selection_pointer != 4))
    {
        cout<<"Invalid entry, please try again: ";
        cin.clear();
        cin.ignore();
        cin>>*year_selection_pointer;

    }
    
    // Displaying selected year:
    // Here i store all the codes of the year my user wants into a vector called auxiliar_year (store index)
    if(*year_selection_pointer != 0)
    {
      for (int j{0}; j < all_courses_and_data.size(); j++)
      {
        if ((int)(std::to_string(course_code[j])[0]) - 48 == *year_selection_pointer)
        {
          auxiliar_year.push_back(j);
        }
      }
      
      //We construc the sub vectors
      for(int a{0}; a < auxiliar_year.size(); a++)
      {
          sub_vector_all_courses_and_data.push_back(all_courses_and_data[auxiliar_year[a]]);
          sub_courses_codes.push_back(course_code[auxiliar_year[a]]);
          sub_courses_names.push_back(course_names[auxiliar_year[a]]);
          sub_grades.push_back(grades[auxiliar_year[a]]);
      }
      // Calculations for selected data
      *mean_of_grades_pointer = mean_calculation(sub_grades);
      *std_deviation_pointer = standard_deviation(sub_grades);
      *error_std_dev_pointer = standard_deviation_error(*std_deviation_pointer, sub_grades.size());

      // Here we display de selected year courses and codes
      cout<<"Sort alphabetically or by code? (a / c): ";
      cin>>*alphabet_or_code_pointer;
      while(cin.fail() or (*alphabet_or_code_pointer != 'a' and *alphabet_or_code_pointer != 'c'))
      {
        cout<<"Invalid entry, please try again: ";
        cin.clear();
        cin.ignore();
        cin>>*alphabet_or_code_pointer;
      }
      
      if(*alphabet_or_code_pointer == 'c')
      {

        std::sort(sub_courses_codes.begin(), sub_courses_codes.end());
        for(int c{0}; c < sub_courses_codes.size(); c++)
        {
          
          for(int d{0}; d < sub_courses_codes.size(); d++)
          {

            if(sub_vector_all_courses_and_data[c].find(std::to_string(sub_courses_codes[d])) != string::npos)
            {
              indexes_code_sort.push_back(d);
            }

          }
        }

        //DISPLAY CODE ORDERED
        for(int f{0}; f < sub_vector_all_courses_and_data.size(); f++)
        {
          for(int g{0}; g < sub_vector_all_courses_and_data.size(); g++)
          {
            if(indexes_code_sort[g] == f)
            {
              ordered_all.push_back(sub_vector_all_courses_and_data[g]);
            }
          }
        }
      // From here alphabetic
      } else if(*alphabet_or_code_pointer == 'a')
      {
        
        std::sort(sub_courses_names.begin(), sub_courses_names.end());
        for(int c{0}; c < sub_courses_names.size(); c++)
        {
          for(int d{0}; d < sub_courses_names.size(); d++)
          {
            if(sub_vector_all_courses_and_data[c].find(sub_courses_names[d]) != string::npos)
            {
              indexes_name_sort.push_back(d);
              break;
            }
          }
        }

        //DISPLAY NAME ORDERED
        for(int f{0}; f < sub_vector_all_courses_and_data.size(); f++)
        {
          for(int g{0}; g < sub_vector_all_courses_and_data.size(); g++)
          {
            if(indexes_name_sort[g] == f)
            {
              ordered_all.push_back(sub_vector_all_courses_and_data[g]);
            }
          }
        }
      } 
    }

    //THE USER CHOSE 0 AND WE SHOULD WORK WITH ALL COURSES
    if(*year_selection_pointer == 0)
    {
      // Making copies of the full course
      course_code_copy = course_code;
      grades_copy = grades;
      course_names_copy = course_names;
      all_courses_and_data_copy = all_courses_and_data;

      cout<<"Sort alphabetically or by code? (a / c): ";
      cin>>*alphabet_or_code_pointer;
      while (cin.fail() or (*alphabet_or_code_pointer != 'a' and *alphabet_or_code_pointer != 'c'))
      {
        cout<<"Invalid entry, please try again: ";
        cin.clear();
        cin.ignore();
        cin>>*alphabet_or_code_pointer;
      }

      if(*alphabet_or_code_pointer == 'c') 
      {
        std::sort(course_code_copy.begin(), course_code_copy.end());
        for(int c{0}; c < course_code_copy.size(); c++)
        {
          for(int d{0}; d < course_code_copy.size(); d++)
          {
            if(all_courses_and_data_copy[c].find(std::to_string(course_code_copy[d])) != string::npos)
            {
              indexes_code_sort.push_back(d);
            }
          }
        }
        //DISPLAY CODE ORDER
        for (int f{0}; f < all_courses_and_data_copy.size(); f++)
        {
          for(int g{0}; g < all_courses_and_data_copy.size(); g++)
          {
            if(indexes_code_sort[g] == f)
              {
                ordered_all.push_back(all_courses_and_data_copy[g]);
              }
          }
        }
      
      } else if (*alphabet_or_code_pointer == 'a')
      {
        std::sort(course_names_copy.begin(), course_names_copy.end());
          for(int c{0}; c < course_names_copy.size(); c++)
          {
            for(int d{0}; d < course_names_copy.size(); d++)
            {
              if(all_courses_and_data_copy[c].find(course_names_copy[d]) != string::npos)
              {
                indexes_name_sort.push_back(d);
                break;
              }
            }
          }

          //DISPLAY NAME ORDERED
          for(int f{0}; f < all_courses_and_data_copy.size(); f++)
          {
            for(int g{0}; g < all_courses_and_data_copy.size(); g++)
            {
              if(indexes_name_sort[g] == f)
              {
                ordered_all.push_back(all_courses_and_data_copy[g]);
              }
            }
          }
      }
    }

    // Printing the results
    cout<<"---"<<endl;
    for(int i{0}; i < ordered_all.size(); i++)
    {
      cout<<ordered_all[i]<<endl;
    }
    cout<<"---"<<endl;
    cout<<"The mean of the grades for the selection is: "<<*mean_of_grades_pointer<<endl;
    cout<<"The standard deviation is: "<<*std_deviation_pointer<<endl;
    cout<<"The standard error in the mean is: "<<*error_std_dev_pointer<<endl;
    cout<<"---"<<endl;

    // Continue running?
    cout<<"Would you like to continue? (y/n): ";
    cin>>*pointer_yes_no;
    while(cin.fail() or (*pointer_yes_no != 'n' and *pointer_yes_no != 'y'))
    {
      cout<<"Invalid entry, please try again: ";
      cin.clear();
      cin.ignore();
      cin>>*pointer_yes_no;
    }

    // Free memory for next run
    sub_vector_all_courses_and_data.clear();
    sub_courses_codes.clear();
    sub_courses_names.clear();
    sub_grades.clear();
    auxiliar_year.clear();
    indexes_code_sort.clear();
    indexes_name_sort.clear();
    ordered_all.clear();
  }

  // Free memory before terminating
  data_file.clear();
  all_courses_and_data.clear();
  grades.clear();
  course_code.clear();
  course_names.clear();
  auxiliar.clear();
  auxiliar_year.clear();
  sub_vector_all_courses_and_data.clear();
  sub_courses_codes.clear();
  sub_courses_names.clear();
  sub_grades.clear();
  indexes_code_sort.clear();
  indexes_name_sort.clear();
  ordered_all.clear();
  all_courses_and_data_copy.clear();
  course_code_copy.clear();
  grades_copy.clear();
  course_names_copy.clear();
  delete point_cont;
  delete course_number_pointer;
  delete std_deviation_pointer;
  delete mean_of_grades_pointer;
  delete error_std_dev_pointer;
  delete pointer_yes_no;
  delete year_selection_pointer;
  delete alphabet_or_code_pointer;

  return 0;
}