#include <argparse.h>

void get_opts(int argc,
              char **argv,
              struct options_t *opts)
{
    if (argc == 1)
    {
        std::cout << "Usage:" << std::endl;
        std::cout << "\t--in or -i <file_path>" << std::endl;
        std::cout << "\t--out or -o <file_path>" << std::endl;
        std::cout << "\t--steps or -s <num_steps>" << std::endl;
        std::cout << "\t--theta or -t <threshold>" << std::endl;
        std::cout << "\t--timestep or -d <timestep>" << std::endl;
        std::cout << "\t[Optional] --visial or -v" << std::endl;
        std::cout << "\t[Optional] --optimize or -p" << std::endl;
        exit(0);
    }

    opts->visual = false;
    opts->op = false;

    struct option l_opts[] = {
        {"in", required_argument, NULL, 'i'},
        {"out", required_argument, NULL, 'o'},
        {"n_steps", required_argument, NULL, 's'},
        {"theta", required_argument, NULL, 't'},
        {"timestep", required_argument, NULL, 'd'},
        {"visial", no_argument, NULL, 'v'},
        {"optimize", no_argument, NULL, 'p'}
    };

    int ind, c;
    while ((c = getopt_long(argc, argv, "i:o:s:t:d:vp", l_opts, &ind)) != -1)
    {
        switch (c)
        {
        case 0:
            break;
        case 'i':
            opts->in_file = (char *)optarg;
            break;
        case 'o':
            opts->out_file = (char *)optarg;
            break;
        case 's':
            opts->n_steps = atoi((char *)optarg);
            break;
        case 't':
            opts->theta = atof((char *)optarg);
            break;
        case 'd':
            opts->time_step = atof((char *)optarg);
            break;
        case 'v':
            opts->visual = true;
            break;
        case 'p':
            opts->op = true;
            break;
        case ':':
            std::cerr << argv[0] << ": option -" << (char)optopt << "requires an argument." << std::endl;
            exit(1);
        }
    }
    
}