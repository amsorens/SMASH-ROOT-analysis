/*
 *  Copyright (c) 2026
 *  Agnieszka Sorensen
 */

//////////////////////////////////////////////////////////////////////////////////////////
// This file contains functions used to plot histograms and Tgraphs
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMASH_ROOT_ANALYSIS_PLOTS_H
#define SMASH_ROOT_ANALYSIS_PLOTS_H

#include <string>
#include <locale>
#include <sstream>

#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TCanvas.h>



std::string double_to_string_with_fixed_decimal_places
  (double value, int number_of_decimal_places);

std::string double_to_string_with_comma(double value, int number_of_decimal_places);

std::string double_to_string_with_the_word_point
  (double value, int number_of_decimal_places);

const char* strip_output_prefix(const char* input);

void set_canvas_properties(TCanvas* can,
			   bool setLogx = false, bool setLogy = false,
			   double left_margin = 0.15, double right_margin = 0.05,
			   double top_margin = 0.08, double bottom_margin = 0.125);

void set_virtual_pad_properties(TVirtualPad* pad,
				int show_ticks_x = 1, int show_ticks_y = 1,
				double left_margin = 0.15, double right_margin = 0.08,
				double top_margin = 0.08, double bottom_margin = 0.15);

void set_axis_properties(TAxis* axis, const char* title, double min = 0.0,
			 double max = 0.0);

void set_2D_histogram_properties(TH2D* histogram_2D,
				 const char* x_axis_label = "",
				 const char* y_axis_label = "");

// Create .C, .pdf, .png files
void create_canvas_files(TCanvas* can,
			 char* basic_file_name);

void create_a_histogram_ROOT_file(TH1D histogram,
				  char* basic_file_name);

// overload of the above for TH1D*
void create_a_histogram_ROOT_file(TH1D* histogram,
				  char* basic_file_name);

// overload of the above for std::unique_ptr<TH1D> (passed by reference to avoid copy)
void create_a_histogram_ROOT_file(const std::unique_ptr<TH1D>& histogram,
				  char* basic_file_name);

// overload of the above for TProfile*
void create_a_histogram_ROOT_file(TProfile* histogram,
				  char* basic_file_name);

// overload of the above for std::unique_ptr<TProfile> (passed by reference to avoid copy)
void create_a_histogram_ROOT_file(const std::unique_ptr<TProfile>& histogram,
				  char* basic_file_name);

// overload of the above for TH2D
void create_a_histogram_ROOT_file(TH2D histogram,
				  char* basic_file_name);

// overload of the above for TH2D*
void create_a_histogram_ROOT_file(TH2D* histogram,
				  char* basic_file_name);

// overload of the above for TH3D*
void create_a_histogram_ROOT_file(TH3D* histogram,
				  char* basic_file_name);

void create_a_tgraph_ROOT_file(TGraphErrors* tgraphs,
			       char* basic_file_name);

void create_a_tgraph_ROOT_file(std::vector< std::vector<TGraphErrors*> > tgraphs,
			       char* basic_file_name);

// overload of the above for multiple tgraphs that are not pointers
void create_a_tgraph_ROOT_file(std::vector<TGraphErrors> tgraphs,
			       char* basic_file_name);

// overload of the above for multiple tgraphs
void create_a_tgraph_ROOT_file(std::vector<TGraphErrors*> tgraphs,
			       char* basic_file_name);

void plot_and_save_1D_histogram(TH1D histogram_1D,
				char* basic_file_name,
				double histogram_range_lower,
				double histogram_range_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy);

void plot_and_save_1D_histogram (TH1D* histogram_1D,
				 char* basic_file_name,
				 double histogram_range_lower,
				 double histogram_range_upper,
				 const char* x_axis_label,
				 const char* y_axis_label = "",
				 const char* plot_option = "h",
				 bool bool_setLogy = "false");

void plot_and_save_1D_histogram(const std::unique_ptr<TH1D>& histogram_1D,
				char* basic_file_name,
				double histogram_range_lower,
				double histogram_range_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy);

void plot_and_save_2D_histogram(TH2D histogram_2D,
				char* basic_file_name,
				double histogram_range_x_lower,
				double histogram_range_x_upper,
				double histogram_range_y_lower,
				double histogram_range_y_upper,
				const char* x_axis_label,
				const char* y_axis_label,
				const char* plot_option,
				bool bool_setLogy);

void plot_and_save_2D_histogram (TH2D* histogram_2D,
				 char* basic_file_name,
				 double histogram_range_x_lower,
				 double histogram_range_x_upper,
				 double histogram_range_y_lower,
				 double histogram_range_y_upper,
				 const char* x_axis_label = "",
				 const char* y_axis_label = "",
				 const char* plot_option = "scat",
				 bool bool_setLogy = "false");


void plot_and_save_3D_histogram (TH3D* histogram_3D,
				 char* basic_file_name,
				 double histogram_range_x_lower,
				 double histogram_range_x_upper,
				 double histogram_range_y_lower,
				 double histogram_range_y_upper,
				 double histogram_range_z_lower,
				 double histogram_range_z_upper,
				 const char* x_axis_label = "",
				 const char* y_axis_label = "",
				 const char* z_axis_label = "",
				 const char* plot_option = "scat",
				 bool bool_setLogy = "false");


void plot_and_save_TGraph(TGraphErrors tgraph,
			  char* basic_file_name,
			  const char* x_axis_label = "",
			  const char* y_axis_label = "",
			  bool bool_setLogx = false,
			  bool bool_setLogy = false,
			  bool use_provided_axis_ranges = false,
			  double x_min_value = -100.0, double x_max_value = -100.0,
			  double y_min_value = -100.0, double y_max_value = -100.0);


void plot_and_save_TGraph(std::vector<TGraphErrors*> tgraph,
			  char* basic_file_name,
			  const char* x_axis_label = "",
			  const char* y_axis_label = "",
			  std::vector<const char*> legend_labels = {},
			  bool bool_setLogx = false,
			  bool bool_setLogy = false,
			  bool use_provided_axis_ranges = false,
			  double x_min_value = -100.0, double x_max_value = -100.0,
			  double y_min_value = -100.0, double y_max_value = -100.0);


#endif
