std::string find_SBMLs (void)
{
    int number_of_SBMLs = 2
    std::vector<std::string> SBML_list[number_of_SBMLs] = {"CAF_Toy_Model.xml","CRC_Toy_Model.xml"}
    for (int i = 0; i < 2; i++)
        std::cout << SBML_list[i] << "\n";
    return SBML_list;
}



void Parse_SBML(void)
{
    
}

