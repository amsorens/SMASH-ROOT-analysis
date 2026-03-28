/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains functions used to plot histograms and Tgraphs
//////////////////////////////////////////////////////////////////////////////////////////

#include "./plots.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>

#include <TFile.h>
#include <TLegend.h>

#include "constants.h"



std::string double_to_string_with_fixed_decimal_places
  (double value, int number_of_decimal_places) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(number_of_decimal_places) << value;
  return oss.str();
}



std::string double_to_string_with_comma(double value, int number_of_decimal_places) {
  // Get the integer and fractional part of a double
  double integer = 0.0;
  double fractional = 0.0;
  fractional = modf(value, &integer);
  const int int_integer = static_cast<int>(integer);
  const int int_fractional =
    static_cast<int>(std::pow(10, number_of_decimal_places) * fractional);
  // Figure out how many leading zeros there are in the fractional
  std::string zeros = "";
  if (fractional > 0 ) {
    double comparison = 0.1;
    for (int i = 0; i < number_of_decimal_places; i++) {
      comparison = std::pow(comparison, i+1);
      if (fractional < comparison) {
	zeros += '0';
      } else {
	break;
      }
    }
  }
  
  char comma_separated[256];
  if (int_fractional == 0) {
    snprintf(comma_separated, sizeof(comma_separated), "%d",	   
	     int_integer);  
  } else {
    snprintf(comma_separated, sizeof(comma_separated), "%d,%s%d",	   
	     int_integer, zeros.c_str(), int_fractional);
  }
  std::string string_output(comma_separated);
  return string_output;
}



std::string double_to_string_with_the_word_point
  (double value, int number_of_decimal_places) {
  // Check the sign
  bool is_negative = value < 0;
  value = std::abs(value);

  // Round to the target number of decimal places
  const double scale = std::pow(10.0, number_of_decimal_places);
  value = std::round(value * scale) / scale;
  
  // Get the integer and fractional part of the double
  double integer = 0.0;
  double fractional = 0.0;
  fractional = std::modf(value, &integer);
  const int int_integer = static_cast<int>(integer);
  const int int_fractional =
    static_cast<int>(std::round(fractional * scale));

  // Format fractional part with leading zeros if needed
  char fractional_with_leading_zeros[Char_Array_Size];
  if (number_of_decimal_places > 0) {
    snprintf(fractional_with_leading_zeros, Char_Array_Size,
	     // this will pad with zeros on the left if needed
	     "%0*d", number_of_decimal_places, int_fractional);
  } else {
    fractional_with_leading_zeros[0] = '\0';
  }
  
  char comma_separated[Char_Array_Size];
  if (number_of_decimal_places == 0) {
    snprintf(comma_separated, Char_Array_Size, "%d", int_integer);  
  } else {
    snprintf(comma_separated, Char_Array_Size, "%dpoint%s",	   
	     int_integer, fractional_with_leading_zeros);
  }
  
  std::string string_output(comma_separated);
  if (is_negative) {
    string_output = "minus" + string_output;
  }
  return string_output;
}



const char* strip_output_prefix(const char* input) {
  // static so that the pointer remains valid after return; note that it will be
  // ovewritten at each use
  static std::string result;
  const std::string prefix = "../output/";

  std::string s(input);

  if (s.find(prefix) == 0) {
    result = s.substr(prefix.length());
  } else {
    result = s;
  }

  return result.c_str();
}



void set_canvas_properties(TCanvas* can,
			   bool setLogx,
			   bool setLogy,
			   double left_margin, double right_margin,
			   double top_margin, double bottom_margin) {
  can->Range(0,0,1,1);
  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetBorderSize(2);
  can->SetLeftMargin(left_margin);
  can->SetRightMargin(right_margin);
  can->SetTopMargin(top_margin);
  can->SetBottomMargin(bottom_margin);
  can->SetFrameBorderMode(0);
  can->SetFrameBorderMode(0);
  if ( setLogx ) {
    can->SetLogx();
  }
  if ( setLogy ) {
    can->SetLogy();
  }
}



void set_virtual_pad_properties(TVirtualPad* pad,
				int show_ticks_x, int show_ticks_y,
				double left_margin, double right_margin,
				double top_margin, double bottom_margin) {
  pad->SetTickx(show_ticks_x);
  pad->SetTicky(show_ticks_y);
  pad->SetLeftMargin(left_margin);
  pad->SetRightMargin(right_margin);
  pad->SetTopMargin(top_margin);
  pad->SetBottomMargin(bottom_margin);
}



void set_axis_properties(TAxis* axis, const char* title, double min, double max) {
  axis->SetTitle(title);
  axis->CenterTitle();
  axis->SetTitleSize(0.055);
  axis->SetTitleOffset(1.4);
  // Set custom limits only if provided by the user
  if ( min != max ) {
    axis->SetLimits(min, max);
  } 
}



void set_2D_histogram_properties(TH2D* histogram_2D,
				 const char* x_axis_label,
				 const char* y_axis_label) {

  histogram_2D->SetStats(0);
  // x axis
  histogram_2D->GetXaxis()->SetTitle(x_axis_label);
  histogram_2D->GetXaxis()->CenterTitle(true);
  //histogram_2D->GetXaxis()->SetTitleSize(0.05);
  histogram_2D->GetXaxis()->SetTitleFont(42); // good options: 42, 132
  histogram_2D->GetXaxis()->SetTitleOffset(1.2);
  histogram_2D->GetXaxis()->SetLabelFont(42);
  //histogram_2D->GetYaxis()->SetLabelSize(0.05);
  // y axis
  histogram_2D->GetYaxis()->SetTitle(y_axis_label);
  //histogram_2D->GetYaxis()->CenterTitle(true);
  //histogram_2D->GetYaxis()->SetTitleSize(0.05);
  histogram_2D->GetYaxis()->SetTitleFont(42); // good options: 42, 132
  histogram_2D->GetYaxis()->SetTitleOffset(1.2);
  histogram_2D->GetYaxis()->SetLabelFont(42);
  //histogram_2D->GetYaxis()->SetLabelSize(0.05); 
}



// Create .C, .pdf, .png files
void create_canvas_files(TCanvas* can,
			 char* basic_file_name) {
  const int char_array_size = 256;

  char histogram_file_name_pdf[char_array_size];
  char histogram_file_name_png[char_array_size];
  char histogram_file_name_C[char_array_size];
  snprintf(histogram_file_name_pdf, char_array_size, "%s.pdf", basic_file_name);
  snprintf(histogram_file_name_png, char_array_size, "%s.png", basic_file_name);
  snprintf(histogram_file_name_C, char_array_size, "%s.C", basic_file_name);

  can->Print(histogram_file_name_pdf);
  can->Print(histogram_file_name_png);
  can->Print(histogram_file_name_C); 
}



void create_a_histogram_ROOT_file(TH1D histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram.Write( histogram.GetName() );

  // close the file
  file->Close();  
}



// overload of the above for TH1D*
void create_a_histogram_ROOT_file(TH1D* histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram->Write( histogram->GetName() );

  // close the file
  file->Close();  
}



// overload of the above for std::unique_ptr<TH1D> (passed by reference to avoid copy)
void create_a_histogram_ROOT_file(const std::unique_ptr<TH1D>& histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram->Write( histogram->GetName() );

  // close the file
  file->Close();  
}



// overload of the above for TProfile*
void create_a_histogram_ROOT_file(TProfile* histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram->Write( histogram->GetName() );

  // close the file
  file->Close();  
}


// overload of the above for std::unique_ptr<TProfile> (passed by reference to avoid copy)
void create_a_histogram_ROOT_file(const std::unique_ptr<TProfile>& histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram->Write( histogram->GetName() );

  // close the file
  file->Close();  
}



// overload of the above for TH2D
void create_a_histogram_ROOT_file(TH2D histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram.Write( histogram.GetName() );

  // close the file
  file->Close();  
}



// overload of the above for TH2D*
void create_a_histogram_ROOT_file(TH2D* histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram->Write( histogram->GetName() );

  // close the file
  file->Close();  
}



// overload of the above for TH3D
void create_a_histogram_ROOT_file(TH3D* histogram,
				  char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char histogram_file_name_root[char_array_size];
  snprintf(histogram_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (histogram_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  histogram->Write( histogram->GetName() );

  // close the file
  file->Close();  
}



void create_a_tgraph_ROOT_file(TGraphErrors* tgraph,
			       char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char tgraph_file_name_root[char_array_size];
  snprintf(tgraph_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (tgraph_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  tgraph->Write( tgraph->GetTitle() );

  // close the file
  file->Close();  
}



// overload of the above for multiple tgraphs
void create_a_tgraph_ROOT_file(std::vector<TGraphErrors*> tgraphs,
			       char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char tgraph_file_name_root[char_array_size];
  snprintf(tgraph_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (tgraph_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  const int number_of_tgraphs = tgraphs.size();
  for (int i = 0; i < number_of_tgraphs; i++) {
    tgraphs[i]->Write( tgraphs[i]->GetTitle() );
  }

  // close the file
  file->Close();  
}



// overload of the above for multiple tgraphs that are not pointers
void create_a_tgraph_ROOT_file(std::vector<TGraphErrors> tgraphs,
			       char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char tgraph_file_name_root[char_array_size];
  snprintf(tgraph_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (tgraph_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  const int number_of_tgraphs = tgraphs.size();
  for (int i = 0; i < number_of_tgraphs; i++) {
    tgraphs[i].Write( tgraphs[i].GetTitle() );
  }

  // close the file
  file->Close();  
}



// overload of the above for multiple multiple tgraphs
void create_a_tgraph_ROOT_file(std::vector< std::vector<TGraphErrors*> > tgraphs,
			       char* basic_file_name) {
  const int char_array_size = 256;

  // Create a ROOT file
  char tgraph_file_name_root[char_array_size];
  snprintf(tgraph_file_name_root, char_array_size, "%s.root", basic_file_name);
  
  TFile *file = new TFile (tgraph_file_name_root, "RECREATE");  
  if ( file->IsOpen() ){
    printf("ROOT file %s opened successfully\n\n", basic_file_name);
  }
  // gFile is a global variable with access to files, which can be accessed by Write()
  gFile = file;

  // objects can be written into the file using their Write method
  const int number_of_tgraphs_1 = tgraphs.size();
  const int number_of_tgraphs_2 = tgraphs[0].size();
  for (int i = 0; i < number_of_tgraphs_1; i++) {
    for (int j = 0; j < number_of_tgraphs_2; j++) {
      tgraphs[i][j]->Write( tgraphs[i][j]->GetTitle() );
    }
  }

  // close the file
  file->Close();  
}



void plot_and_save_1D_histogram(TH1D histogram_1D,
				char* basic_file_name,
				double histogram_range_lower,
				double histogram_range_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy) {
  ///////////////////////////////////////////
  // Plot the histogram
  
  //declare a canvas with name, title, size x, size y
  TCanvas *can = new TCanvas("can","Canvas",1000,1000);
  set_canvas_properties(can, bool_setLogy);

  histogram_1D.GetXaxis()->SetRangeUser(histogram_range_lower,
					histogram_range_upper);
  histogram_1D.SetMarkerColor(kBlack);
  histogram_1D.SetTitle("");

  histogram_1D.GetXaxis()->SetTitle(x_axis_label);
  histogram_1D.GetXaxis()->CenterTitle(true);
  histogram_1D.GetXaxis()->SetTitleSize(0.05);
  histogram_1D.GetXaxis()->SetLabelFont(132);
  histogram_1D.GetXaxis()->SetLabelSize(0.05);
  histogram_1D.GetXaxis()->SetTitleOffset(1.2);
  histogram_1D.GetXaxis()->SetTitleFont(132);

  histogram_1D.GetYaxis()->SetTitle(y_axis_label);
  histogram_1D.GetYaxis()->CenterTitle(true);
  histogram_1D.GetYaxis()->SetTitleSize(0.05);
  histogram_1D.GetYaxis()->SetLabelFont(132);
  histogram_1D.GetYaxis()->SetLabelSize(0.05);
  histogram_1D.GetYaxis()->SetTitleOffset(1.2);
  histogram_1D.GetYaxis()->SetTitleFont(132);
  
  histogram_1D.Draw(plot_option);


  
  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(can, basic_file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_histogram_ROOT_file(histogram_1D, basic_file_name);
}



void plot_and_save_1D_histogram(TH1D* histogram_1D,
				char* basic_file_name,
				double histogram_range_lower,
				double histogram_range_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy) {
  ///////////////////////////////////////////
  // Plot the histogram
  
  //declare a canvas with name, title, size x, size y
  TCanvas *can = new TCanvas("can","Canvas",1000,1000);
  set_canvas_properties(can, bool_setLogy);

  histogram_1D->GetXaxis()->SetRangeUser(histogram_range_lower,
					 histogram_range_upper);
  histogram_1D->SetMarkerColor(kBlack);
  histogram_1D->SetTitle("");

  histogram_1D->GetXaxis()->SetTitle(x_axis_label);
  histogram_1D->GetXaxis()->CenterTitle(true);
  histogram_1D->GetXaxis()->SetTitleSize(0.05);
  histogram_1D->GetXaxis()->SetLabelFont(132);
  histogram_1D->GetXaxis()->SetLabelSize(0.05);
  histogram_1D->GetXaxis()->SetTitleOffset(1.2);
  histogram_1D->GetXaxis()->SetTitleFont(132);

  histogram_1D->GetYaxis()->SetTitle(y_axis_label);
  histogram_1D->GetYaxis()->CenterTitle(true);
  histogram_1D->GetYaxis()->SetTitleSize(0.05);
  histogram_1D->GetYaxis()->SetLabelFont(132);
  histogram_1D->GetYaxis()->SetLabelSize(0.05);
  histogram_1D->GetYaxis()->SetTitleOffset(1.2);
  histogram_1D->GetYaxis()->SetTitleFont(132);
  
  histogram_1D->Draw(plot_option);


  
  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(can, basic_file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_histogram_ROOT_file(histogram_1D, basic_file_name);
}



void plot_and_save_1D_histogram(const std::unique_ptr<TH1D>& histogram_1D,
				char* basic_file_name,
				double histogram_range_lower,
				double histogram_range_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy) {
  ///////////////////////////////////////////
  // Plot the histogram
  
  //declare a canvas with name, title, size x, size y
  TCanvas *can = new TCanvas("can","Canvas",1000,1000);
  set_canvas_properties(can, bool_setLogy);

  histogram_1D->GetXaxis()->SetRangeUser(histogram_range_lower,
					 histogram_range_upper);
  histogram_1D->SetMarkerColor(kBlack);
  histogram_1D->SetTitle("");

  histogram_1D->GetXaxis()->SetTitle(x_axis_label);
  histogram_1D->GetXaxis()->CenterTitle(true);
  histogram_1D->GetXaxis()->SetTitleSize(0.05);
  histogram_1D->GetXaxis()->SetLabelFont(132);
  histogram_1D->GetXaxis()->SetLabelSize(0.05);
  histogram_1D->GetXaxis()->SetTitleOffset(1.2);
  histogram_1D->GetXaxis()->SetTitleFont(132);

  histogram_1D->GetYaxis()->SetTitle(y_axis_label);
  histogram_1D->GetYaxis()->CenterTitle(true);
  histogram_1D->GetYaxis()->SetTitleSize(0.05);
  histogram_1D->GetYaxis()->SetLabelFont(132);
  histogram_1D->GetYaxis()->SetLabelSize(0.05);
  histogram_1D->GetYaxis()->SetTitleOffset(1.2);
  histogram_1D->GetYaxis()->SetTitleFont(132);
  
  histogram_1D->Draw(plot_option);


  
  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(can, basic_file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_histogram_ROOT_file(histogram_1D, basic_file_name);
}



void plot_and_save_2D_histogram(TH2D histogram_2D,
				char* basic_file_name,
				double histogram_range_x_lower,
				double histogram_range_x_upper,
				double histogram_range_y_lower,
				double histogram_range_y_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy) {
  ////////////////////////////
  // Plot the histogram
  
  //declare a canvas with name, title, size x, size y
  TCanvas *can = new TCanvas("can","Canvas",1000,1000);

  set_canvas_properties(can, bool_setLogy);

  histogram_2D.GetXaxis()->SetRangeUser(histogram_range_x_lower,
					histogram_range_x_upper);
  histogram_2D.GetYaxis()->SetRangeUser(histogram_range_y_lower,
					histogram_range_y_upper);

  histogram_2D.SetMarkerColor(kBlack);
  histogram_2D.SetTitle("");

  histogram_2D.GetXaxis()->SetTitle(x_axis_label);
  histogram_2D.GetXaxis()->CenterTitle(true);
  histogram_2D.GetXaxis()->SetTitleSize(0.05);
  histogram_2D.GetXaxis()->SetLabelFont(132);
  histogram_2D.GetXaxis()->SetLabelSize(0.05);
  histogram_2D.GetXaxis()->SetTitleOffset(1.2);
  histogram_2D.GetXaxis()->SetTitleFont(132);

  histogram_2D.GetYaxis()->SetTitle(y_axis_label);
  histogram_2D.GetYaxis()->CenterTitle(true);
  histogram_2D.GetYaxis()->SetTitleSize(0.05);
  histogram_2D.GetYaxis()->SetLabelFont(132);
  histogram_2D.GetYaxis()->SetLabelSize(0.05);
  histogram_2D.GetYaxis()->SetTitleOffset(1.2);
  histogram_2D.GetYaxis()->SetTitleFont(132);
  
  histogram_2D.Draw(plot_option);


  
  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(can, basic_file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_histogram_ROOT_file(histogram_2D, basic_file_name);
}



// overload of the above for a TH2D*
void plot_and_save_2D_histogram(TH2D* histogram_2D,
				char* basic_file_name,
				double histogram_range_x_lower,
				double histogram_range_x_upper,
				double histogram_range_y_lower,
				double histogram_range_y_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy) {
  ////////////////////////////
  // Plot the histogram
  
  //declare a canvas with name, title, size x, size y
  TCanvas *can = new TCanvas("can","Canvas",1000,1000);

  set_canvas_properties(can, bool_setLogy);

  histogram_2D->GetXaxis()->SetRangeUser(histogram_range_x_lower,
					 histogram_range_x_upper);
  histogram_2D->GetYaxis()->SetRangeUser(histogram_range_y_lower,
					 histogram_range_y_upper);

  histogram_2D->SetMarkerColor(kBlack);
  histogram_2D->SetTitle("");

  histogram_2D->GetXaxis()->SetTitle(x_axis_label);
  histogram_2D->GetXaxis()->CenterTitle(true);
  histogram_2D->GetXaxis()->SetTitleSize(0.05);
  histogram_2D->GetXaxis()->SetLabelFont(132);
  histogram_2D->GetXaxis()->SetLabelSize(0.05);
  histogram_2D->GetXaxis()->SetTitleOffset(1.2);
  histogram_2D->GetXaxis()->SetTitleFont(132);

  histogram_2D->GetYaxis()->SetTitle(y_axis_label);
  histogram_2D->GetYaxis()->CenterTitle(true);
  histogram_2D->GetYaxis()->SetTitleSize(0.05);
  histogram_2D->GetYaxis()->SetLabelFont(132);
  histogram_2D->GetYaxis()->SetLabelSize(0.05);
  histogram_2D->GetYaxis()->SetTitleOffset(1.2);
  histogram_2D->GetYaxis()->SetTitleFont(132);
  
  histogram_2D->Draw(plot_option);


  
  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(can, basic_file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_histogram_ROOT_file(histogram_2D, basic_file_name);
}



void plot_and_save_3D_histogram(TH3D* histogram_3D,
				char* basic_file_name,
				double histogram_range_x_lower,
				double histogram_range_x_upper,
				double histogram_range_y_lower,
				double histogram_range_y_upper,
				double histogram_range_z_lower,
				double histogram_range_z_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* z_axis_label,
				const char* plot_option,
				bool bool_setLogy) {
  ////////////////////////////
  // Plot the histogram
  
  //declare a canvas with name, title, size x, size y
  TCanvas *can = new TCanvas("can","Canvas",1000,1000);

  set_canvas_properties(can, bool_setLogy);

  histogram_3D->GetXaxis()->SetRangeUser(histogram_range_x_lower,
					 histogram_range_x_upper);
  histogram_3D->GetYaxis()->SetRangeUser(histogram_range_y_lower,
					 histogram_range_y_upper);
  histogram_3D->GetZaxis()->SetRangeUser(histogram_range_z_lower,
					 histogram_range_z_upper);

  histogram_3D->SetMarkerColor(kBlack);
  histogram_3D->SetTitle("");

  histogram_3D->GetXaxis()->SetTitle(x_axis_label);
  histogram_3D->GetXaxis()->CenterTitle(true);
  histogram_3D->GetXaxis()->SetTitleSize(0.05);
  histogram_3D->GetXaxis()->SetLabelFont(132);
  histogram_3D->GetXaxis()->SetLabelSize(0.05);
  histogram_3D->GetXaxis()->SetTitleOffset(1.2);
  histogram_3D->GetXaxis()->SetTitleFont(132);

  histogram_3D->GetYaxis()->SetTitle(y_axis_label);
  histogram_3D->GetYaxis()->CenterTitle(true);
  histogram_3D->GetYaxis()->SetTitleSize(0.05);
  histogram_3D->GetYaxis()->SetLabelFont(132);
  histogram_3D->GetYaxis()->SetLabelSize(0.05);
  histogram_3D->GetYaxis()->SetTitleOffset(1.2);
  histogram_3D->GetYaxis()->SetTitleFont(132);

  histogram_3D->GetZaxis()->SetTitle(z_axis_label);
  histogram_3D->GetZaxis()->CenterTitle(true);
  histogram_3D->GetZaxis()->SetTitleSize(0.05);
  histogram_3D->GetZaxis()->SetLabelFont(132);
  histogram_3D->GetZaxis()->SetLabelSize(0.05);
  histogram_3D->GetZaxis()->SetTitleOffset(1.2);
  histogram_3D->GetZaxis()->SetTitleFont(132);
  
  histogram_3D->Draw(plot_option);


  
  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(can, basic_file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_histogram_ROOT_file(histogram_3D, basic_file_name);
}



void plot_and_save_TGraph(TGraphErrors tgraph,
			  char* basic_file_name,
			  const char* x_axis_label,
			  const char* y_axis_label,
			  bool bool_setLogx,
			  bool bool_setLogy,
			  bool use_provided_axis_ranges,
			  double x_min_value, double x_max_value,
			  double y_min_value, double y_max_value) {
  ///////////////////////////////////////////
  // Plot the TGraph

  ///////////////////////////////////////////
  // Declare a canvas with name, title, size x, size y
  std::unique_ptr<TCanvas> tgraph_canvas
    = std::make_unique<TCanvas>(basic_file_name, basic_file_name,
				200, 200, 500, 500);
  set_canvas_properties(tgraph_canvas.get(), bool_setLogx, bool_setLogy);
  
  ///////////////////////////////////////////
  // Establish min, max
  
  // Initialize variables that will store the minimum and maximum values to absurd values
  double x_min = 1000000000;
  double y_min = 1000000000;
  double x_max = -1000000000;
  double y_max = -1000000000;
  double x_error_max = -1000000000;
  double y_error_max = -1000000000;

  // Use provided values
  if ( use_provided_axis_ranges ) {
    x_min = x_min_value;
    x_max = x_max_value;
    y_min = y_min_value;
    y_max = y_max_value;
  }
  // Figure out the values
  else {
    const int n_points = tgraph.GetN();
    // Loop over all points in the TGraph
    for (int j = 0; j < n_points; j++) {
      // Get x and y coordinates of the point
      const double x = tgraph.GetX()[j];
      const double y = tgraph.GetY()[j];
      const double x_error = tgraph.GetEX()[j];
      const double y_error = tgraph.GetEY()[j];
      // Update minimum and maximum values as appropriate
      if (x < x_min) { x_min = x; }
      if (x > x_max) { x_max = x; }
      if (y < y_min) { y_min = y; }
      if (y > y_max) { y_max = y; }
      if (x_error > x_error_max) { x_error_max = x_error; }
      if (y_error > y_error_max) { y_error_max = y_error; }
    }

    if (!bool_setLogx && !bool_setLogy) {
      // Redefine the max values to include the errors
      x_min = x_min - x_error_max;
      x_max = x_max + x_error_max;
      y_min = y_min - y_error_max;
      y_max = y_max + y_error_max;
    }

    // Define the spread in values
    const double x_spread = x_max - x_min;
    const double y_spread = y_max - y_min;

    x_min = bool_setLogx ? (x_min / 5.0) : (x_min - 0.2 * x_spread);
    x_max = bool_setLogx ? (x_max * 5.0) : (x_max + 0.2 * x_spread);
    y_min = bool_setLogy ? (y_min / 5.0) : (y_min - 0.2 * y_spread);
    y_max = bool_setLogy ? (y_max * 5.0) : (y_max + 0.2 * y_spread);
  }
  
  

  ///////////////////////////////////////////
  // Define an auxiliary histogram to control the plotting area
  char auxiliary_histogram_name[Char_Array_Size];
  snprintf(auxiliary_histogram_name, Char_Array_Size,
	   "%s_auxiliary_h", tgraph.GetName());
  // Clean the basic file name for the histogram title
  const char* auxiliary_histogram_title = strip_output_prefix(basic_file_name);
  std::unique_ptr<TH2D> h_auxiliary = std::make_unique<TH2D>
    (auxiliary_histogram_name, auxiliary_histogram_title,
     100, x_min, x_max, 100, y_min, y_max);
  // Set properties for a nice plot
  set_2D_histogram_properties(h_auxiliary.get(), x_axis_label, y_axis_label);

  ///////////////////////////////////////////
  // Draw the histogram, the TGraph, and the fit
  h_auxiliary->Draw();
  //tgraph.SetLineColor(0);
  tgraph.SetMarkerStyle(21);
  //tgraph.SetMarkerColor(kRed);
  tgraph.Draw("psame");
  tgraph_canvas->Print();

  ///////////////////////////////////////////
  // Update the name if needed
  char file_name[Char_Array_Size];
  if (bool_setLogx || bool_setLogy) { 
    if (bool_setLogx && bool_setLogy) {
      snprintf(file_name, Char_Array_Size,
	       "%s_log_log", basic_file_name);
    }
    else {
      snprintf(file_name, Char_Array_Size,
	       "%s_log", basic_file_name);
    }
  } else {
    snprintf(file_name, Char_Array_Size,
	     "%s", basic_file_name);
  }

  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(tgraph_canvas.get(), file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_tgraph_ROOT_file(&tgraph, file_name);
  
  ///////////////////////////////////////////
  // Delete pointers?
  //delete chi2_canvas;
}



void plot_and_save_TGraph(std::vector<TGraphErrors*> tgraph,
			  char* basic_file_name,
			  const char* x_axis_label,
			  const char* y_axis_label,
			  std::vector<const char*> legend_labels,
			  bool bool_setLogx,
			  bool bool_setLogy,
			  bool use_provided_axis_ranges,
			  double x_min_value, double x_max_value,
			  double y_min_value, double y_max_value) {
  ///////////////////////////////////////////
  // Plot the TGraph

  ///////////////////////////////////////////
  // Declare a canvas with name, title, size x, size y
  std::unique_ptr<TCanvas> tgraph_canvas
    = std::make_unique<TCanvas>(basic_file_name, basic_file_name,
				200, 200, 500, 500);
  set_canvas_properties(tgraph_canvas.get(), bool_setLogx, bool_setLogy);
  
  ///////////////////////////////////////////
  // Establish min, max
  
  // Initialize variables that will store the minimum and maximum values to absurd values
  double x_min = 1000000000;
  double y_min = 1000000000;
  double x_max = -1000000000;
  double y_max = -1000000000;
  double x_error_max = -1000000000;
  double y_error_max = -1000000000;

  // Use provided values
  if ( use_provided_axis_ranges ) {
    x_min = x_min_value;
    x_max = x_max_value;
    y_min = y_min_value;
    y_max = y_max_value;
  }
  // Figure out the values
  else {
    // Loop over all TGraphs
    for (int i = 0; i < static_cast<int>(tgraph.size()); i++) { 
      const int n_points = tgraph[i]->GetN();
      // Loop over all points in the TGraph
      for (int j = 0; j < n_points; j++) {
	// Get x and y coordinates of the point
	const double x = tgraph[i]->GetX()[j];
	const double y = tgraph[i]->GetY()[j];
	const double x_error = tgraph[i]->GetEX()[j];
	const double y_error = tgraph[i]->GetEY()[j];
	// Update minimum and maximum values as appropriate
	if (x < x_min) { x_min = x; }
	if (x > x_max) { x_max = x; }
	if (y < y_min) { y_min = y; }
	if (y > y_max) { y_max = y; }
	if (x_error > x_error_max) { x_error_max = x_error; }
	if (y_error > y_error_max) { y_error_max = y_error; }
      }
    }

    if (!bool_setLogx && !bool_setLogy) {
      // Redefine the max values to include the errors
      x_min = x_min - x_error_max;
      x_max = x_max + x_error_max;
      y_min = y_min - y_error_max;
      y_max = y_max + y_error_max;
    }

    // Define the spread in values
    const double x_spread = x_max - x_min;
    const double y_spread = y_max - y_min;

    x_min = bool_setLogx ? (x_min / 5.0) : (x_min - 0.2 * x_spread);
    x_max = bool_setLogx ? (x_max * 5.0) : (x_max + 0.2 * x_spread);
    y_min = bool_setLogy ? (y_min / 5.0) : (y_min - 0.2 * y_spread);
    y_max = bool_setLogy ? (y_max * 5.0) : (y_max + 0.2 * y_spread);
  }
  
  

  ///////////////////////////////////////////
  // Define an auxiliary histogram to control the plotting area
  char auxiliary_histogram_name[Char_Array_Size];
  snprintf(auxiliary_histogram_name, Char_Array_Size,
	   "%s_auxiliary_h", tgraph[0]->GetName());
  // Clean the basic file name for the histogram title
  const char* auxiliary_histogram_title = strip_output_prefix(basic_file_name);
  std::unique_ptr<TH2D> h_auxiliary = std::make_unique<TH2D>
    (auxiliary_histogram_name, auxiliary_histogram_title,
     100, x_min, x_max, 100, y_min, y_max);
  // Set properties for a nice plot
  set_2D_histogram_properties(h_auxiliary.get(), x_axis_label, y_axis_label);

  ///////////////////////////////////////////
  // Helper vector of colors
  std::vector<Color_t> colors =
    {kRed, kBlue, kBlack, kMagenta, kGreen+3, kOrange+7,kViolet-3};

  ///////////////////////////////////////////
  // Draw the histogram and the TGraph
  h_auxiliary->Draw();
  // Create the legend (adjust coordinates as needed: x1,y1,x2,y2)
  TLegend legend(0.65, 0.7, 0.88, 0.88);
  legend.SetBorderSize(0); // Remove the border
  legend.SetFillStyle(0); // Transparent background
  for (int i = 0; i < static_cast<int>(tgraph.size()); i++) {
    tgraph[i]->SetLineColor(colors[i]);
    tgraph[i]->SetMarkerStyle(21 + i);
    tgraph[i]->SetMarkerColor(colors[i]);
    tgraph[i]->Draw("psame");

    if (i < static_cast<int>(legend_labels.size()) && legend_labels[i]) {
      legend.AddEntry(tgraph[i], legend_labels[i], "lp");
    }
  }
  legend.Draw("same");
  tgraph_canvas->Print();

  ///////////////////////////////////////////
  // Update the name if needed
  char file_name[Char_Array_Size];
  if (bool_setLogx || bool_setLogy) { 
    if (bool_setLogx && bool_setLogy) {
      snprintf(file_name, Char_Array_Size,
	       "%s_log_log", basic_file_name);
    }
    else {
      snprintf(file_name, Char_Array_Size,
	       "%s_log", basic_file_name);
    }
  } else {
    snprintf(file_name, Char_Array_Size,
	     "%s", basic_file_name);
  }

  ///////////////////////////////////////////
  // Create .C, .pdf, .png files
  create_canvas_files(tgraph_canvas.get(), file_name);

  ///////////////////////////////////////////
  // Create a ROOT file
  create_a_tgraph_ROOT_file(tgraph, file_name);
  
  ///////////////////////////////////////////
  // Delete pointers?
  //delete chi2_canvas;
}
