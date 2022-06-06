#include <gtk/gtk.h>
#include "gtkchart.h"

static void print_hello (GtkWidget *widget,gpointer   data){
  g_print ("Hello World\n");
}

static void activate (GtkApplication *app, gpointer user_data){
  GtkWidget *window;
  GtkWidget *button;
  gtk_chart_get_type();

  GtkWidget *chart = gtk_chart_new();
  gtk_chart_set_type((GtkChart*) chart, GTK_CHART_TYPE_SCATTER);
  gtk_chart_set_logscale_y((GtkChart*) chart,10);
  gtk_chart_set_title((GtkChart*)chart, "Title");
  gtk_chart_set_label((GtkChart*)chart, "Label");
  gtk_chart_set_x_label((GtkChart*)chart, "X label [ ]");
  gtk_chart_set_y_label((GtkChart*)chart, "Y label [ ]");
  gtk_chart_set_x_max((GtkChart*)chart, 3);
  gtk_chart_set_y_max((GtkChart*)chart, 1000);
  gtk_chart_set_width((GtkChart*)chart, 800);
  gtk_chart_plot_point((GtkChart*)chart, 0.0, 0.0);
  gtk_chart_plot_point((GtkChart*)chart, 1.0, 10.0);
  gtk_chart_plot_point((GtkChart*)chart, 2.0, 50.0);
  gtk_chart_plot_point((GtkChart*)chart, 3.0, 1000.0);

  window = gtk_application_window_new (app);
  gtk_window_set_title (GTK_WINDOW (window), "Window");
  gtk_window_set_default_size (GTK_WINDOW (window), 200, 200);

  button = gtk_button_new_with_label ("Hello World");
  g_signal_connect (button, "clicked", G_CALLBACK (print_hello), NULL);
 // gtk_window_set_child (GTK_WINDOW (window), button);
  gtk_window_set_child (GTK_WINDOW (window), chart);

  gtk_window_present (GTK_WINDOW (window));
}

int main (int argc,char **argv){
  GtkApplication *app;
  int status;

  app = gtk_application_new ("org.gtk.example", G_APPLICATION_FLAGS_NONE);
  g_signal_connect (app, "activate", G_CALLBACK (activate), NULL);
  status = g_application_run (G_APPLICATION (app), argc, argv);
  g_object_unref (app);

  return status;
}
