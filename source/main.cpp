// PROJECT includes
#include <wave/parameters/ParameterHandler.hh>
#include <wave/ElasticWave_cGcG.tpl.hh>
#include <wave/ElasticWave_cG.tpl.hh>
#include <wave/ElasticWave_dG.tpl.hh>
#include <wave/Problem/Problem.hh>


// DEAL.II includes



// C++ includes
#include <iostream>
#include <fstream>
#include <memory>
#include <filesystem>
#include <sys/stat.h>  // for mkdir
#include <limits>


namespace fs = std::filesystem;

void clear_output_directory(const std::string& output_path) {
    if (fs::exists(output_path)) {
        for (const auto& entry : fs::directory_iterator(output_path)) {
            fs::remove_all(entry.path());
        }
    }
}


void display_menu(const std::vector<std::string>& options, const std::string& prompt) {
    int space = prompt.length();
    std::cout << std::string(76, '=') << "\n";
    std::cout << "=" << std::string(74, ' ') << "=\n";
    std::cout << "=" << std::string((74-space)/2, ' ') << prompt << std::string((74-space)/2, ' ') << "=\n";
    std::cout << "=" << std::string(74, ' ') << "=\n";
    std::cout << std::string(76, '=') << "\n";
    for (size_t i = 0; i < options.size(); ++i) {
        std::cout << "[ " << i + 1 << " ] " << options[i] << "\n";
    }
    std::cout << std::string(76, '=') << "\n";
    std::cout << "Your choice: ";
}

size_t get_user_choice(size_t max_choice) {
    size_t choice;
    while (true) {
        std::cin >> choice;
        if (std::cin.fail() || choice < 1 || choice > max_choice) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid selection, try again: ";
        } else {
            break;
        }
    }
    return choice;
}

std::string choose_folder(const std::vector<std::string>& folders) {
    display_menu(folders, "Choose a method:");
    return folders[get_user_choice(folders.size()) - 1];
}

std::string choose_file(const std::string& folder) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.path().extension() == ".prm") {
            files.push_back(entry.path().filename().string());
        }
    }
    
    if (files.empty()) {
        std::cerr << "Did not find any .prm files." << std::endl;
        exit(1);
    }
    std::sort(files.begin(), files.end());
    display_menu(files, "Choose one example: ");
    return folder + "/" + files[get_user_choice(files.size()) - 1];
}




int main() {

    try
    {
        dealii::deallog.depth_console(2);

        std::string input_folder = "../input"; // Basisordner
        std::vector<std::string> subfolders = {"cG", "dG", "cG-cG"};
        
        std::string selected_folder = choose_folder(subfolders);
        std::string full_path = input_folder + "/" + selected_folder;
        
        if (!fs::exists(full_path)) {
            std::cerr << "Method does not exsist: " << full_path << std::endl;
            return 1;
        }
        
        std::string selected_file = choose_file(full_path);
        
        std::cout << "\nUsing input file: " << selected_file << "\n";

        std::cout << "============================================================================" << std::endl;
        std::cout << "=                                                                          =" << std::endl;
        std::cout << "=                               RUNNING EXAMPLE                            =" << std::endl;
        std::cout << "=                                                                          =" << std::endl;
        std::cout << "============================================================================" << std::endl;

        auto parameter_handler = std::make_shared<wave::ParameterHandler>();
        parameter_handler->parse_input(selected_file);

        int dim = parameter_handler->get_integer("dim");
        unsigned int s = parameter_handler->get_integer("s");
        unsigned int r = parameter_handler->get_integer("r");
        
        std::string output_path = "../output/" + selected_folder + "/" + fs::path(selected_file).stem().string() + "/"
                                  + "final_condition=" + std::to_string(parameter_handler->get_bool("Junker Boundary"));

        clear_output_directory(output_path);
        
        std::shared_ptr < wave::Problem > problem;

        
        switch (dim) {
            case 1: {
                if (selected_folder == "cG-cG")
                    problem = std::make_shared< ElasticWave_cGcG<1> > (s, r);
                else if (selected_folder == "dG")
                    problem = std::make_shared< ElasticWave_dG<1> > (s, r);
                else if (selected_folder == "cG")
                    problem = std::make_shared< ElasticWave_cG<1> > (s, r);
                // SpaceTimeFEM<1> problem;
                break;
            }
            case 2: {
                if (selected_folder == "cG-cG")
                    problem = std::make_shared< ElasticWave_cGcG<2> > (s, r);
                else if (selected_folder == "dG")
                    problem = std::make_shared< ElasticWave_dG<2> > (s, r);
                else if (selected_folder == "cG")
                    problem = std::make_shared< ElasticWave_cG<2> > (s, r);
                // SpaceTimeFEM<1> problem;
                break;
            }

            default: {
                dealii::ExcNotImplemented();
                break;
            }
        }

        problem->set_input_parameters(parameter_handler);
        problem->run();

    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;

        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }

    return 0;
}