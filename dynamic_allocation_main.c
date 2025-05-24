#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>

// edge id'ler ve demand id'ler 1den başlıyor

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define MAX_LINE_LENGTH 256
#define MUTATION_RATE 0.15 // %10 mutasyon oranı
#define NUMBER_OF_NODES 1243
#define NUMBER_OF_EDGES 9138
#define MAX_TIME_STRING_LENGTH 9


#define NUMBER_OF_LANES 2

// store alternative path's for the OD's
#define MAX_NUMBER_OF_EDGES_IN_ROUTE 100
#define NUMBER_OF_ALTERNATIVES 5
#define NUMBER_OF_ODS 11 * 14

// store the edge matrix
#define NUMBER_OF_TIME_INTERVALS 1200
#define NUMBER_OF_LINES 2
#define POPULATION_SIZE 50 // Population sizd must be even
#define MULTIPLICITY_STRING_REP 10
#define NUMBER_OF_GENERATIONS 20000

// Katsayılar sonrasında öneme göre değişecek
#define W1 0.5
#define W2 0.5
#define W3 0.5
#define PENALTY 10

#define CROSSOVER_TYPE 3 // 0: order, 1: problem specific, 2: similarity crossover minimize conflict, 3: similarity crossover fitness related
#define MUTATION_TYPE 4  // 0: type1, 1: type2, 2: type3, 3: type4, 4: type5

typedef struct uatfm_route_network_edge
{
    int edge_id;
    int start_node;
    int end_node;
    int start_layer;
    int end_layer;
    int edge_type;
    float length;
    float total_cost;
    float flight_cost;
    float noise_cost;
} UATFM_Route_Network_Edge;

typedef struct string_rep
{
    short od_sequence[NUMBER_OF_ODS * MULTIPLICITY_STRING_REP];
    char selected_alternative[NUMBER_OF_ODS * MULTIPLICITY_STRING_REP];
    int max_delays[NUMBER_OF_ODS * MULTIPLICITY_STRING_REP];
    int real_delays[NUMBER_OF_ODS * MULTIPLICITY_STRING_REP];
} String_Rep;

typedef struct node
{
    int node_id;
    int layer; // On which layer nodes is located
} Node;

typedef struct od_demand
{
    Node source_node;      // Source node of the OD pair
    Node destination_node; // Destination node of the OD pair
    char *EST;             // Earliest Start Time
    char *LAT;             // Latest Arrival Time
    int UAV_type;          // UAV type to indicate velocity
    int start_edge_id;
    int end_edge_id;
    int priority; // number 1: means highest priority & number 3 : means lowest priority
} OD_Demand;

OD_Demand *demands;

int number_of_demand;

UATFM_Route_Network_Edge *UATFM_route_network;
unsigned short ***all_alternative_paths_for_ODs;
unsigned short ***alternative_paths_for_ODs;
int ****edge_matrix;
int ***edge_matrix_for_individual1;
int ***edge_matrix_for_individual2;
int *edges_for_f3;
String_Rep *population;
String_Rep *population_backup;

double compute_fitness(String_Rep solution, int solution_index);
void evaluate_population_fitness(String_Rep *population, double *fitness_values);
int tournament_selection(String_Rep *population, double *fitness_values);
void crossover(String_Rep parent1, String_Rep parent2, String_Rep *child1, String_Rep *child2);
void shuffle(short *array, int n);
OD_Demand *read_od_demand(const char *filename);
int is_blank_line(const char *line);
void type1_mutation(String_Rep *population);
double compute_demand_fitness_contribution(String_Rep solution, int od_index);
void fill_edge_matrix_for_individual(String_Rep *rep, int individual_index);
void low_similarity_crossover(String_Rep parent1, String_Rep parent2, String_Rep *child1, String_Rep *child2, String_Rep *population, int parent1_index, int parent2_index);

void allocate_all()
{
    UATFM_route_network = malloc(sizeof(UATFM_Route_Network_Edge) * (NUMBER_OF_EDGES + 1));
    if (!UATFM_route_network)
    {
        fprintf(stderr, "[ERROR] malloc failed for UATFM_route_network\n");
        exit(EXIT_FAILURE);
    }
    // Burada eğer tüm struct sıfırlanacaksa:
    memset(UATFM_route_network, 0, sizeof(UATFM_Route_Network_Edge) * (NUMBER_OF_EDGES + 1));

    edges_for_f3 = malloc(sizeof(int) * NUMBER_OF_EDGES);
    if (!edges_for_f3)
    {
        fprintf(stderr, "[ERROR] malloc failed for edges_for_f3\n");
        exit(EXIT_FAILURE);
    }
    memset(edges_for_f3, 0, sizeof(int) * NUMBER_OF_EDGES);

    population = calloc(POPULATION_SIZE, sizeof(String_Rep));
    if (!population)
    {
        fprintf(stderr, "[ERROR] calloc failed for population\n");
        exit(EXIT_FAILURE);
    }

    population_backup = calloc(POPULATION_SIZE, sizeof(String_Rep));
    if (!population_backup)
    {
        fprintf(stderr, "[ERROR] calloc failed for population_backup\n");
        exit(EXIT_FAILURE);
    }

    all_alternative_paths_for_ODs = malloc(NUMBER_OF_ODS * sizeof(unsigned short **));
    if (!all_alternative_paths_for_ODs)
    {
        fprintf(stderr, "[ERROR] malloc failed for all_alternative_paths_for_ODs\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < NUMBER_OF_ODS; i++)
    {
        all_alternative_paths_for_ODs[i] = malloc((NUMBER_OF_ALTERNATIVES + 1) * sizeof(unsigned short *));
        if (!all_alternative_paths_for_ODs[i])
        {
            fprintf(stderr, "[ERROR] malloc failed for all_alternative_paths_for_ODs[%d]\n", i);
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j <= NUMBER_OF_ALTERNATIVES; j++)
        {
            all_alternative_paths_for_ODs[i][j] = calloc(MAX_NUMBER_OF_EDGES_IN_ROUTE, sizeof(unsigned short));
            if (!all_alternative_paths_for_ODs[i][j])
            {
                fprintf(stderr, "[ERROR] calloc failed for all_alternative_paths_for_ODs[%d][%d]\n", i, j);
                exit(EXIT_FAILURE);
            }
            // calloc zaten sıfırlıyor, ekstra memset gerek yok
        }
    }

    alternative_paths_for_ODs = malloc(NUMBER_OF_ODS * MULTIPLICITY_STRING_REP * sizeof(unsigned short **));
    if (!alternative_paths_for_ODs)
    {
        fprintf(stderr, "[ERROR] malloc failed for alternative_paths_for_ODs\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < NUMBER_OF_ODS * MULTIPLICITY_STRING_REP; i++)
    {
        alternative_paths_for_ODs[i] = malloc(NUMBER_OF_ALTERNATIVES * sizeof(unsigned short *));
        if (!alternative_paths_for_ODs[i])
        {
            fprintf(stderr, "[ERROR] malloc failed for alternative_paths_for_ODs[%d]\n", i);
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < NUMBER_OF_ALTERNATIVES; j++)
        {
            alternative_paths_for_ODs[i][j] = calloc(MAX_NUMBER_OF_EDGES_IN_ROUTE, sizeof(unsigned short));
            if (!alternative_paths_for_ODs[i][j])
            {
                fprintf(stderr, "[ERROR] calloc failed for alternative_paths_for_ODs[%d][%d]\n", i, j);
                exit(EXIT_FAILURE);
            }
            // calloc zaten sıfırlıyor, ekstra memset gerek yok
        }
    }

    edge_matrix = malloc(POPULATION_SIZE * sizeof(int ***));
    if (!edge_matrix)
    {
        fprintf(stderr, "[ERROR] malloc failed for edge_matrix\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        edge_matrix[i] = malloc(NUMBER_OF_TIME_INTERVALS * sizeof(int **));
        if (!edge_matrix[i])
        {
            fprintf(stderr, "[ERROR] malloc failed for edge_matrix[%d]\n", i);
            exit(EXIT_FAILURE);
        }
        for (int t = 0; t < NUMBER_OF_TIME_INTERVALS; t++)
        {
            edge_matrix[i][t] = malloc(NUMBER_OF_EDGES * sizeof(int *));
            if (!edge_matrix[i][t])
            {
                fprintf(stderr, "[ERROR] malloc failed for edge_matrix[%d][%d]\n", i, t);
                exit(EXIT_FAILURE);
            }
            for (int e = 0; e < NUMBER_OF_EDGES; e++)
            {
                edge_matrix[i][t][e] = malloc(NUMBER_OF_LINES * sizeof(int));
                if (!edge_matrix[i][t][e])
                {
                    fprintf(stderr, "[ERROR] malloc failed for edge_matrix[%d][%d][%d]\n", i, t, e);
                    exit(EXIT_FAILURE);
                }
                // memset ile hızlı -1 ataması
                memset(edge_matrix[i][t][e], 0xFF, NUMBER_OF_LINES * sizeof(int));
            }
        }
    }

    edge_matrix_for_individual1 = malloc(NUMBER_OF_TIME_INTERVALS * sizeof(int **));
    edge_matrix_for_individual2 = malloc(NUMBER_OF_TIME_INTERVALS * sizeof(int **));
    if (!edge_matrix_for_individual1 || !edge_matrix_for_individual2)
    {
        fprintf(stderr, "[ERROR] malloc failed for edge_matrix_for_individual1/2\n");
        exit(EXIT_FAILURE);
    }
    for (int t = 0; t < NUMBER_OF_TIME_INTERVALS; t++)
    {
        edge_matrix_for_individual1[t] = malloc(NUMBER_OF_EDGES * sizeof(int *));
        edge_matrix_for_individual2[t] = malloc(NUMBER_OF_EDGES * sizeof(int *));
        if (!edge_matrix_for_individual1[t] || !edge_matrix_for_individual2[t])
        {
            fprintf(stderr, "[ERROR] malloc failed for edge_matrix_for_individual1/2[%d]\n", t);
            exit(EXIT_FAILURE);
        }
        for (int e = 0; e < NUMBER_OF_EDGES; e++)
        {
            edge_matrix_for_individual1[t][e] = calloc(NUMBER_OF_LINES, sizeof(int));
            edge_matrix_for_individual2[t][e] = calloc(NUMBER_OF_LINES, sizeof(int));
            if (!edge_matrix_for_individual1[t][e] || !edge_matrix_for_individual2[t][e])
            {
                fprintf(stderr, "[ERROR] calloc failed for edge_matrix_for_individual1/2[%d][%d]\n", t, e);
                exit(EXIT_FAILURE);
            }
            // calloc zaten sıfırlıyor, ekstra memset gerek yok
        }
    }
}

void free_all()
{
    free(UATFM_route_network);
    free(edges_for_f3);
    free(population);
    free(population_backup);

    for (int i = 0; i < NUMBER_OF_ODS; i++)
    {
        for (int j = 0; j <= NUMBER_OF_ALTERNATIVES; j++)
        {
            free(all_alternative_paths_for_ODs[i][j]);
        }
        free(all_alternative_paths_for_ODs[i]);
    }
    free(all_alternative_paths_for_ODs);

    for (int i = 0; i < NUMBER_OF_ODS * MULTIPLICITY_STRING_REP; i++)
    {
        for (int j = 0; j < NUMBER_OF_ALTERNATIVES; j++)
        {
            free(alternative_paths_for_ODs[i][j]);
        }
        free(alternative_paths_for_ODs[i]);
    }
    free(alternative_paths_for_ODs);

    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        for (int t = 0; t < NUMBER_OF_TIME_INTERVALS; t++)
        {
            for (int e = 0; e < NUMBER_OF_EDGES; e++)
            {
                free(edge_matrix[i][t][e]);
            }
            free(edge_matrix[i][t]);
        }
        free(edge_matrix[i]);
    }
    free(edge_matrix);

    for (int t = 0; t < NUMBER_OF_TIME_INTERVALS; t++)
    {
        for (int e = 0; e < NUMBER_OF_EDGES; e++)
        {
            free(edge_matrix_for_individual1[t][e]);
            free(edge_matrix_for_individual2[t][e]);
        }
        free(edge_matrix_for_individual1[t]);
        free(edge_matrix_for_individual2[t]);
    }
    free(edge_matrix_for_individual1);
    free(edge_matrix_for_individual2);
}

void subtract_times(char *time1, char *time2, char *result)
{
    int hour1, minute1, second1;
    int hour2, minute2, second2;

    // get hour minute and second using "HH:MM:SS"
    sscanf(time1, "%d:%d:%d", &hour1, &minute1, &second1);
    sscanf(time2, "%d:%d:%d", &hour2, &minute2, &second2);

    // if the first time is greater than the second time, add 24 hours to get the correct difference
    int totalSeconds1 = hour1 * 3600 + minute1 * 60 + second1;
    int totalSeconds2 = hour2 * 3600 + minute2 * 60 + second2;

    int diffSeconds = totalSeconds1 - totalSeconds2;

    if (diffSeconds < 0)
    {
        diffSeconds += 24 * 3600; // if negative, add 24 hours
    }

    // convert the difference into hours, minutes and seconds
    int diffHours = diffSeconds / 3600;
    int diffMinutes = (diffSeconds % 3600) / 60;
    int diffSecs = diffSeconds % 60;

    // format the result into a string
    sprintf(result, "%02d:%02d:%02d", diffHours, diffMinutes, diffSecs);
}

int get_string_length(String_Rep rep)
{
    int length = 0;
    while (rep.od_sequence[length] != 0)
    {
        length++;
    }
    return length;
}

void add_seconds_to_time(const char *original, int seconds_to_add, char *result)
{
    int h, m, s;
    sscanf(original, "%d:%d:%d", &h, &m, &s);

    int total_seconds = h * 3600 + m * 60 + s + seconds_to_add;

    h = (total_seconds / 3600) % 24;
    m = (total_seconds % 3600) / 60;
    s = total_seconds % 60;

    sprintf(result, "%02d:%02d:%02d", h, m, s);
}

void subtract_hours_from_time(char *time, double hoursToSubtract, char *result)
{
    int hour, minute, second;

    // get the hour minute and second from "HH:MM:SS"
    sscanf(time, "%d:%d:%d", &hour, &minute, &second);

    // convert it into seconds
    int total_seconds = hour * 3600 + minute * 60 + second;

    // convert the time to subtract into seconds
    int subtractSeconds = (int)(hoursToSubtract * 3600);

    // subtract from the total seconds
    total_seconds -= subtractSeconds;

    // to handle negqtive values, add 24 hours and take the modulo
    total_seconds = (total_seconds % 86400 + 86400) % 86400;

    // convert it into hours, minutes and seconds
    hour = total_seconds / 3600;
    minute = (total_seconds % 3600) / 60;
    second = total_seconds % 60;

    // format it into a string and write it to result
    sprintf(result, "%02d:%02d:%02d", hour, minute, second);
}

int string_time_to_seconds(const char *timeDiff)
{
    int hours, minutes, seconds;

    // get hours, minutes and seconds from "HH:MM:SS"
    sscanf(timeDiff, "%d:%d:%d", &hours, &minutes, &seconds);

    // calculate total seconds
    return (hours * 3600) + (minutes * 60) + seconds;
}

int find_max_delay(String_Rep solution, int index)
{
    // get the solution's od sequence at the specified index
    int od = solution.od_sequence[index] - 1; // -1 to match the array index
    // get the selected alternative of the solution at the specified index
    int selected = solution.selected_alternative[index];

    int edge;
    double total_length = 0;

    // loop through the edges
    for (int i = 0; i < NUMBER_OF_EDGES; i++)
    {
        // if the edge is 0, break the loop
        if (alternative_paths_for_ODs[od][selected][i] == 0)
        {
            break;
        }
        // get the selected edge and add its length to the total length
        edge = alternative_paths_for_ODs[od][selected][i];
        total_length += UATFM_route_network[edge].length;
    }

    // calculate the flight time
    double flight_time = total_length / demands[od].UAV_type;

    // UAV type value in -> km / h
    // flight time in hours

    // LAT - flight time = LST
    // LST - EST = max delay

    char LST[9]; // "HH:MM:SS"
    // subtract the flight time from the LAT and write it to LST
    subtract_hours_from_time(demands[od].LAT, flight_time, LST);

    char max_delay[9];

    if (string_time_to_seconds(LST) < string_time_to_seconds(demands[od].EST))
    {
        printf("Error: There is an inconsistency between the flight time and the latest arrival time; it is not possible for the UAV to arrive at the desired time.\n\
            Error on demand number %d\n",
               od);
        return -1;
    }

    // subtract the EST from the LST and write it to max_delay
    subtract_times(LST, demands[od].EST, max_delay);

    // then we convert it into seconds
    int max_delay_in_sec = string_time_to_seconds(max_delay);

    // check if the max delay is negative
    if (max_delay_in_sec < 0)
    {
        printf("Error: Negative delay value at index %d\n", index);
    }

    return max_delay_in_sec;
}

// Satırın tamamen boş olup olmadığını kontrol eden fonksiyon
int is_blank_line(const char *line)
{
    return (strlen(line) == strspn(line, " \r\n\t"));
}

void cancel_invalid_demands(String_Rep *child, int max_length)
{
    int length = 0;
    int write_index = 0;

    // uzunluk fonksiyonu çağırılacak
    for (int read_index = 0; read_index < max_length; read_index++)
    {
        if (child->od_sequence[read_index] != 0)
        {
            if (!(child->od_sequence[read_index] == -1))
            {
                child->od_sequence[write_index] = child->od_sequence[read_index];
                child->selected_alternative[write_index] = child->selected_alternative[read_index];
                child->max_delays[write_index] = child->max_delays[read_index];
                child->real_delays[write_index] = child->real_delays[read_index];
                write_index++;
            }
        }
        else
        {
            break;
        }
    }

    length = write_index;

    // uzunluk fonksiyonu çağırılacak
    for (int i = length; i < max_length; i++)
    {
        child->max_delays[i] = 0;
        child->od_sequence[i] = 0;
        child->selected_alternative[i] = 0;
        child->real_delays[i] = 0;
    }
}

void shuffle_priority_aware(short *array, int k)
{
    // Öncelik seviyelerine göre ayrı ayrı index listeleri oluştur
    int count[3] = {0}; // priority 1 -> index 0, priority 2 -> 1, priority 3 -> 2
    short *priority_buckets[3];

    // Geçici bellek ayır
    for (int p = 0; p < 3; p++)
    {
        priority_buckets[p] = malloc(sizeof(short) * k); // maksimum k olabilir
    }

    // priority'ye göre bölüştür
    for (int i = 0; i < k; i++)
    {
        int od_index = array[i];
        int p = demands[od_index].priority - 1; // priority 1 => index 0
        if (p >= 0 && p < 3)
        {
            priority_buckets[p][count[p]++] = od_index;
        }
    }

    // Her öncelik grubunu kendi içinde karıştır
    for (int p = 0; p < 3; p++)
    {
        for (int i = count[p] - 1; i > 0; i--)
        {
            int j = rand() % (i + 1);
            short tmp = priority_buckets[p][i];
            priority_buckets[p][i] = priority_buckets[p][j];
            priority_buckets[p][j] = tmp;
        }
    }

    // array'e sırayla geri yerleştir (öncelikli olan başta olacak şekilde)
    int idx = 0;
    for (int p = 0; p < 3; p++)
    {
        for (int i = 0; i < count[p]; i++)
        {
            array[idx++] = priority_buckets[p][i];
        }
    }

    // Belleği temizle
    for (int p = 0; p < 3; p++)
    {
        free(priority_buckets[p]);
    }
}

// Select a parent using tournament selection
int tournament_selection(String_Rep *population, double *fitness_values)
{
    int tournament_size = 5; // Number of individuals in the tournament
    int best_index = -1;
    double best_fitness = DBL_MAX; // Initialize the best fitness value to the maximum possible value

    for (int i = 0; i < tournament_size; i++)
    {
        int random_index = rand() % POPULATION_SIZE; // Randomly select an individual
        if (best_index == -1 || fitness_values[random_index] < best_fitness)
        {
            best_index = random_index;
            best_fitness = fitness_values[random_index];
        }
    }

    return best_index; // Return the best individual
}

double find_upperbound(int od)
{
    double upperbound = 0;
    //  seçili od için alternatifleri dolaşır ve en uzun olanı bulur
    for (int i = 0; i < NUMBER_OF_ALTERNATIVES; i++)
    {
        int selected_edge;
        double length = 0;
        for (int j = 0; j < NUMBER_OF_EDGES; j++)
        {
            if (alternative_paths_for_ODs[od][i][j] == 0)
            {
                break;
            }
            selected_edge = alternative_paths_for_ODs[od][i][j];
            length += UATFM_route_network[selected_edge].length;
        }
        if (length > upperbound)
        {
            upperbound = length;
        }
    }
    return upperbound;
}

int find_unsatisified_demands(String_Rep solution)
{
    int unsatisfied_count = 0;
    int solution_length = get_string_length(solution);

    for (int i = 0; i < solution_length; i++)
    {
        if (solution.real_delays[i] > solution.max_delays[i])
        {
            unsatisfied_count++;
        }
    }

    return unsatisfied_count;
}

// Function to compute fitness value for a single solution
double compute_fitness(String_Rep solution, int solution_index)
{
    double f1 = 0, f2 = 0, f3 = 0;
    int number_of_paths_involved = 0;
    // Find the number of paths involved
    while (solution.od_sequence[number_of_paths_involved] != 0)
    {
        number_of_paths_involved++;
    }

    int solution_length = get_string_length(solution);

    // Compute f1: Path Length Related
    for (int i = 0; i < solution_length; i++)
    {

        // solution içindeki od_sequence ve selected_alternative değerlerini al
        // path length hesapla ve max length ile oranını al
        int od = solution.od_sequence[i] - 1; // -1 to match the array index
        int selected = solution.selected_alternative[i];
        int selected_edge;
        double length = 0;
        double max_length = find_upperbound(od);
        for (int j = 0; j < NUMBER_OF_EDGES; j++)
        {

            if (alternative_paths_for_ODs[od][selected][j] == 0)
            {
                break;
            }

            selected_edge = alternative_paths_for_ODs[od][selected][j];
            length += UATFM_route_network[selected_edge].length;
        }
        f1 += length / max_length;
    }

    f1 /= number_of_paths_involved;

    // Compute f2: Delay Related
    double total_delay_ratio = 0;

    for (int i = 0; i < solution_length; i++)
    {

        if (solution.max_delays[i] == 0)
        {
            printf("Warning: Division by zero at index %d\n", i);
            continue;
        }
        // solution'ın real_delays ve max_delays değerlerini al
        // real_delays / max_delays oranını al
        total_delay_ratio += (double)solution.real_delays[i] / solution.max_delays[i];
    }

    if (number_of_paths_involved != 0)
    {
        f2 = total_delay_ratio / number_of_paths_involved;
    }
    else
    {
        printf("Error: number_of_paths_involved is zero.\n");
    }

    // f3 = (double)nonzero_count / (zero_count + nonzero_count);
    int nonzero_count = 0;
    int total_count = 0;

    memset(edges_for_f3, 0, NUMBER_OF_EDGES * sizeof(int));

    for (int i = 0; i < solution_length; i++)
    {
        for (int j = 0; j < NUMBER_OF_EDGES; j++)
        {
            int od = solution.od_sequence[i] - 1; // -1 to match the array index
            int selected = solution.selected_alternative[i];

            if (alternative_paths_for_ODs[od][selected][j] == 0)
            {
                break;
            }

            int selected_edge = alternative_paths_for_ODs[od][selected][j];

            if (edges_for_f3[selected_edge] == 0)
            {
                edges_for_f3[selected_edge] = 1;
            }
        }
    }

    for (int j = 0; j < NUMBER_OF_EDGES; j++)
    {
        if (edges_for_f3[j] != 0)
        {
            total_count++;
        }
        else
        {
            continue;
        }
    }

    total_count = total_count * NUMBER_OF_TIME_INTERVALS * NUMBER_OF_LANES;

    for (int t = 0; t < NUMBER_OF_TIME_INTERVALS; t++)
    {
        for (int e = 0; e < NUMBER_OF_EDGES; e++)
        {
            for (int l = 0; l < NUMBER_OF_LINES; l++)
            {
                nonzero_count += (edge_matrix[solution_index][t][e][l] != -1);
            }
        }
    }

    f3 = (double)nonzero_count / total_count;

    int unsatisfied_count = find_unsatisified_demands(solution);
    // Compute overall fitness function
    // değerler ve katsayıları çarp ve topla
    double fitness_value = (W1 * f1) + (W2 * f2) + (W3 * f3) + +unsatisfied_count * PENALTY;

    return fitness_value;
}

// Function to compute fitness value for a single solution
double *compute_fitness_for_best_solution(String_Rep solution, int solution_index)
{
    double f1 = 0, f2 = 0, f3 = 0;
    int number_of_paths_involved = 0;
    int solution_length = get_string_length(solution);
    double *fitness_values = malloc(10 * sizeof(double));
    if (fitness_values == NULL)
    {
        fprintf(stderr, "Bellek ayrilamadi.\n");
        exit(1);
    }
    // Find the number of paths involved
    while (solution.od_sequence[number_of_paths_involved] != 0)
    {
        number_of_paths_involved++;
    }

    // Compute f1: Path Length Related
    for (int i = 0; i < solution_length; i++)
    {

        // solution içindeki od_sequence ve selected_alternative değerlerini al
        // path length hesapla ve max length ile oranını al
        int od = solution.od_sequence[i] - 1; // -1 to match the array index
        int selected = solution.selected_alternative[i];
        int selected_edge;
        double length = 0;
        double max_length = find_upperbound(od);
        for (int j = 0; j < NUMBER_OF_EDGES; j++)
        {

            if (alternative_paths_for_ODs[od][selected][j] == 0)
            {
                break;
            }

            selected_edge = alternative_paths_for_ODs[od][selected][j];
            length += UATFM_route_network[selected_edge].length;
        }
        f1 += length / max_length;
    }
    fitness_values[4] = f1;

    f1 /= number_of_paths_involved;
    fitness_values[5] = (double)number_of_paths_involved;

    // Compute f2: Delay Related
    double total_delay_ratio = 0;

    for (int i = 0; i < solution_length; i++)
    {

        if (solution.max_delays[i] == 0)
        {
            printf("Warning: Division by zero at index %d\n", i);
            continue;
        }
        // solution'ın real_delays ve max_delays değerlerini al
        // real_delays / max_delays oranını al
        total_delay_ratio += (double)solution.real_delays[i] / solution.max_delays[i];
    }
    fitness_values[6] = total_delay_ratio;
    if (number_of_paths_involved != 0)
    {
        f2 = total_delay_ratio / number_of_paths_involved;
    }
    else
    {
        printf("Error: number_of_paths_involved is zero.\n");
    }

    // Compute f3: Congestion Related
    // edge_matrix içindeki 0 olmayan eleman sayısını ve 0 olan eleman sayısını bul
    // 0 olmayan eleman sayısının toplam eleman sayısına oranını al
    int nonzero_count = 0;
    int total_count = 0;

    for (int i = 0; i < solution_length; i++)
    {
        for (int j = 0; j < NUMBER_OF_EDGES; j++)
        {
            int od = solution.od_sequence[i] - 1; // -1 to match the array index
            int selected = solution.selected_alternative[i];

            if (alternative_paths_for_ODs[od][selected][j] == 0)
            {
                break;
            }

            int selected_edge = alternative_paths_for_ODs[od][selected][j];

            if (edges_for_f3[selected_edge] == 0)
            {
                edges_for_f3[selected_edge] = 1;
            }
        }
    }

    for (int j = 0; j < NUMBER_OF_EDGES; j++)
    {
        if (edges_for_f3[j] != 0)
        {
            total_count++;
        }
        else
        {
            continue;
        }
    }

    total_count = total_count * NUMBER_OF_TIME_INTERVALS * NUMBER_OF_LANES;

    for (int t = 0; t < NUMBER_OF_TIME_INTERVALS; t++)
    {
        for (int e = 0; e < NUMBER_OF_EDGES; e++)
        {
            for (int l = 0; l < NUMBER_OF_LINES; l++)
            {
                nonzero_count += (edge_matrix[solution_index][t][e][l] != -1);
            }
        }
    }

    fitness_values[7] = (double)nonzero_count;
    fitness_values[8] = (double)total_count;

    f3 = (double)nonzero_count / total_count;

    // Compute overall fitness function
    // değerler ve katsayıları çarp ve topla

    int unsatisfied_count = find_unsatisified_demands(solution);
    fitness_values[9] = (double)unsatisfied_count;
    double fitness_value = (W1 * f1) + (W2 * f2) + (W3 * f3) + unsatisfied_count * PENALTY;

    // double fitness_values[4] = {f1, f2, f3, fitness_value};

    fitness_values[0] = f1;
    fitness_values[1] = f2;
    fitness_values[2] = f3;
    fitness_values[3] = fitness_value;

    return fitness_values;
}

float find_similarity_ratio(String_Rep parent1, String_Rep parent2)
{
    int length_parent_1 = get_string_length(parent1);
    int length_parent_2 = get_string_length(parent2);
    int min_length = MIN(length_parent_1, length_parent_2);
    int max_length = MAX(length_parent_1, length_parent_2);
    int similarity = 0;

    for (int i = 0; i < min_length; i++)
    {
        if (parent1.od_sequence[i] == parent2.od_sequence[i] && parent1.selected_alternative[i] == parent2.selected_alternative[i])
        {
            similarity++;
        }
    }
    return ((float)similarity / max_length) * 100.0f;
}

String_Rep *sort_by_priority(String_Rep *individual)
{
    int length = get_string_length(*individual);
    // Sort the individual based on the priority of the demands
    for (int i = 0; i < length - 1; i++)
    {
        for (int j = 0; j < length - i - 1; j++)
        {
            if (demands[individual->od_sequence[j] - 1].priority > demands[individual->od_sequence[j + 1] - 1].priority)
            {
                // od_sequence'deki elemanları takas et
                short temp = individual->od_sequence[j];
                individual->od_sequence[j] = individual->od_sequence[j + 1];
                individual->od_sequence[j + 1] = temp;

                // selected_alternative'ı takas et
                char temp_alt = individual->selected_alternative[j];
                individual->selected_alternative[j] = individual->selected_alternative[j + 1];
                individual->selected_alternative[j + 1] = temp_alt;

                // max_delays'ı takas et
                int temp_max_delay = individual->max_delays[j];
                individual->max_delays[j] = individual->max_delays[j + 1];
                individual->max_delays[j + 1] = temp_max_delay;

                // real_delays'ı takas et
                int temp_real_delay = individual->real_delays[j];
                individual->real_delays[j] = individual->real_delays[j + 1];
                individual->real_delays[j + 1] = temp_real_delay;
            }
        }
    }
    return individual;
}

void print_string_rep(String_Rep rep)
{
    printf("od_sequence:           ");
    for (int i = 0; i < get_string_length(rep); i++)
        printf("%d ", rep.od_sequence[i]);
    printf("\nselected_alternative:  ");
    for (int i = 0; i < get_string_length(rep); i++)
        printf("%d ", rep.selected_alternative[i]);
    printf("\nmax_delays:            ");
    for (int i = 0; i < get_string_length(rep); i++)
        printf("%d ", rep.max_delays[i]);
    printf("\nreal_delays:           ");
    for (int i = 0; i < get_string_length(rep); i++)
        printf("%d ", rep.real_delays[i]);
    printf("\n\n");
}

void order_crossover(String_Rep parent1, String_Rep parent2, String_Rep *child1, String_Rep *child2)
{
    int start, end, i, j, k;

    int length1 = get_string_length(parent1);
    int length2 = get_string_length(parent2);
    int child_length = MAX(length1, length2);

    // child1 ve child2'yi sıfırla
    for (i = 0; i <= child_length; i++)
    {
        child1->od_sequence[i] = 0;
        child2->od_sequence[i] = 0;
    }

    // Rastgele start ve end seç
    start = rand() % length1;
    end = start + (rand() % (length1 - start));

    // Segmenti kopyala
    for (i = start; i <= end; i++)
    {
        child1->od_sequence[i] = parent1.od_sequence[i];
        child1->selected_alternative[i] = parent1.selected_alternative[i];
        child1->max_delays[i] = parent1.max_delays[i];
        child1->real_delays[i] = parent1.real_delays[i];

        child2->od_sequence[i] = parent2.od_sequence[i];
        child2->selected_alternative[i] = parent2.selected_alternative[i];
        child2->max_delays[i] = parent2.max_delays[i];
        child2->real_delays[i] = parent2.real_delays[i];
    }

    // child1 doldur
    k = (end + 1) % child_length;
    for (i = 0; i < length2; i++)
    {
        int node = parent2.od_sequence[i], exists = 0;
        for (j = start; j <= end; j++)
        {
            if (child1->od_sequence[j] == node)
            {
                exists = 1;
                break;
            }
        }
        if (!exists)
        {
            int loop_counter = 0;
            while (child1->od_sequence[k] != 0 && loop_counter < child_length)
            {
                k = (k + 1) % child_length;
                loop_counter++;
            }
            if (loop_counter == child_length)
                continue;
            child1->od_sequence[k] = node;
            child1->selected_alternative[k] = parent2.selected_alternative[i];
            child1->max_delays[k] = parent2.max_delays[i];
            child1->real_delays[k] = parent2.real_delays[i];
        }
    }
    // parent1'in kalanları da child1'e doldur
    for (i = 0; i < length1; i++)
    {
        int node = parent1.od_sequence[i], exists = 0;
        for (j = 0; j < child_length; j++)
        {
            if (child1->od_sequence[j] == node)
            {
                exists = 1;
                break;
            }
        }
        if (!exists)
        {
            int loop_counter = 0;
            while (child1->od_sequence[k] != 0 && loop_counter < child_length)
            {
                k = (k + 1) % child_length;
                loop_counter++;
            }
            if (loop_counter == child_length)
                continue;
            child1->od_sequence[k] = node;
            child1->selected_alternative[k] = parent1.selected_alternative[i];
            child1->max_delays[k] = parent1.max_delays[i];
            child1->real_delays[k] = parent1.real_delays[i];
        }
    }

    // child2 doldur
    k = (end + 1) % child_length;
    for (i = 0; i < length1; i++)
    {
        int node = parent1.od_sequence[i], exists = 0;
        for (j = start; j <= end; j++)
        {
            if (child2->od_sequence[j] == node)
            {
                exists = 1;
                break;
            }
        }
        if (!exists)
        {
            int loop_counter = 0;
            while (child2->od_sequence[k] != 0 && loop_counter < child_length)
            {
                k = (k + 1) % child_length;
                loop_counter++;
            }
            if (loop_counter == child_length)
                continue;
            child2->od_sequence[k] = node;
            child2->selected_alternative[k] = parent1.selected_alternative[i];
            child2->max_delays[k] = parent1.max_delays[i];
            child2->real_delays[k] = parent1.real_delays[i];
        }
    }
    // parent2'nin kalanlarını child2'ye doldur
    for (i = 0; i < length2; i++)
    {
        int node = parent2.od_sequence[i], exists = 0;
        for (j = 0; j < child_length; j++)
        {
            if (child2->od_sequence[j] == node)
            {
                exists = 1;
                break;
            }
        }
        if (!exists)
        {
            int loop_counter = 0;
            while (child2->od_sequence[k] != 0 && loop_counter < child_length)
            {
                k = (k + 1) % child_length;
                loop_counter++;
            }
            if (loop_counter == child_length)
                continue;
            child2->od_sequence[k] = node;
            child2->selected_alternative[k] = parent2.selected_alternative[i];
            child2->max_delays[k] = parent2.max_delays[i];
            child2->real_delays[k] = parent2.real_delays[i];
        }
    }

    // son 0 koy
    child1->od_sequence[child_length] = 0;
    child2->od_sequence[child_length] = 0;
    child1 = sort_by_priority(child1);
    child2 = sort_by_priority(child2);
}

void similarity_crossover_fitness_related(String_Rep parent1, String_Rep parent2, String_Rep *child)
{
    int max_length = MAX(get_string_length(parent1), get_string_length(parent2));
    int child_index = 0;

    for (int i = 0; i < max_length; i++)
    {
        short od1 = parent1.od_sequence[i];
        short od2 = parent2.od_sequence[i];

        // Parent1 veya Parent2'den biri bitmişse
        if (od1 == 0 || od2 == 0)
        {
            String_Rep *active_parent = (od1 != 0) ? &parent1 : &parent2;

            double contribution = compute_demand_fitness_contribution(*active_parent, i);
            if (contribution <= 5.0)
            {
                child->od_sequence[child_index] = active_parent->od_sequence[i];
                child->selected_alternative[child_index] = active_parent->selected_alternative[i];
                child->max_delays[child_index] = active_parent->max_delays[i];
                child_index++;
            }
            continue;
        }

        // Aynı OD ve alternatifse direkt yaz
        if (od1 == od2 && parent1.selected_alternative[i] == parent2.selected_alternative[i])
        {
            child->od_sequence[child_index] = od1;
            child->selected_alternative[child_index] = parent1.selected_alternative[i];
            child->max_delays[child_index] = parent1.max_delays[i];
            child_index++;
        }
        else
        {
            // Contribution karşılaştırması
            double c1 = compute_demand_fitness_contribution(parent1, i);
            double c2 = compute_demand_fitness_contribution(parent2, i);

            String_Rep *chosen = (c1 <= c2) ? &parent1 : &parent2;

            child->od_sequence[child_index] = chosen->od_sequence[i];
            child->selected_alternative[child_index] = chosen->selected_alternative[i];
            child->max_delays[child_index] = chosen->max_delays[i];
            child_index++;
        }
    }

    // Child dizisini sonlandır
    if (child_index < NUMBER_OF_ODS * MULTIPLICITY_STRING_REP)
    {
        child->od_sequence[child_index] = 0;
    }

    child = sort_by_priority(child);
}

void similarity_crossover_minimize_conflict(String_Rep parent1, String_Rep parent2, String_Rep *child)
{
    int max_length = MAX(get_string_length(parent1), get_string_length(parent2));

    int child_index = 0;
    for (int i = 0; i < max_length; i++)
    {
        if (parent1.od_sequence[i] == 0 || parent2.od_sequence[i] == 0)
        {
            if (parent1.od_sequence[i] != 0)
            {
                if (parent1.real_delays[i] < parent1.max_delays[i])
                {
                    child->od_sequence[child_index] = parent1.od_sequence[i];
                    child->selected_alternative[child_index] = parent1.selected_alternative[i];
                    child->max_delays[child_index] = parent1.max_delays[i];
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if (parent2.real_delays[i] < parent2.max_delays[i])
                {
                    child->od_sequence[child_index] = parent2.od_sequence[i];
                    child->selected_alternative[child_index] = parent2.selected_alternative[i];
                    child->max_delays[child_index] = parent2.max_delays[i];
                }
                else
                {
                    continue;
                }
            }
        }
        // eğer real delayleri max delayi aşıyorsa continue
        else if (parent1.real_delays[i] > parent1.max_delays[i] && parent2.real_delays[i] > parent2.max_delays[i])
        {
            continue;
        }
        // eğer demandlar ve selected alternatif sequence'leri aynı ise
        else if (parent1.od_sequence[i] == parent2.od_sequence[i] && parent1.selected_alternative[i] == parent2.selected_alternative[i])
        {
            child->od_sequence[child_index] = parent1.od_sequence[i];
            child->selected_alternative[child_index] = parent1.selected_alternative[i];
            child->max_delays[child_index] = parent1.max_delays[i];
        }
        // eğer parent1 ve parent2'nin od_sequence'leri aynı değilse ama delaylerinde sıkıntı yoksa
        else if (parent1.real_delays[i] < parent1.max_delays[i] && parent2.real_delays[i] < parent2.max_delays[i])
        {
            // Prioritye göre seç
            if (demands[parent1.od_sequence[i]].priority == demands[parent2.od_sequence[i]].priority)
            {
                if (parent1.real_delays[i] < parent2.real_delays[i])
                {
                    child->od_sequence[child_index] = parent1.od_sequence[i];
                    child->selected_alternative[child_index] = parent1.selected_alternative[i];
                    child->max_delays[child_index] = parent1.max_delays[i];
                }
                else
                {
                    child->od_sequence[child_index] = parent2.od_sequence[i];
                    child->selected_alternative[child_index] = parent2.selected_alternative[i];
                    child->max_delays[child_index] = parent2.max_delays[i];
                }
            }
            else
            {
                // farklı prioritylerdeyle, prioritye göre seç
                if (demands[parent1.od_sequence[i]].priority < demands[parent2.od_sequence[i]].priority)
                {
                    child->od_sequence[child_index] = parent1.od_sequence[i];
                    child->selected_alternative[child_index] = parent1.selected_alternative[i];
                    child->max_delays[child_index] = parent1.max_delays[i];
                }
                else
                {
                    child->od_sequence[child_index] = parent2.od_sequence[i];
                    child->selected_alternative[child_index] = parent2.selected_alternative[i];
                    child->max_delays[child_index] = parent2.max_delays[i];
                }
            }
        }
        // eğer ikisinden birinin real delayi max'ı geçmiyorsa
        else
        {
            if (parent1.real_delays[i] < parent1.max_delays[i])
            {
                child->od_sequence[i] = parent1.od_sequence[i];
                child->selected_alternative[i] = parent1.selected_alternative[i];
                child->max_delays[i] = parent1.max_delays[i];
            }
            else
            {
                child->od_sequence[i] = parent2.od_sequence[i];
                child->selected_alternative[i] = parent2.selected_alternative[i];
                child->max_delays[i] = parent2.max_delays[i];
            }
        }
        child_index++;
    }
}

OD_Demand *read_od_demand(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Dosya açılamadı");
        return NULL;
    }

    OD_Demand *demands = NULL;
    int count = 0; // Toplam OD_Demand sayısı
    int start_id, start_layer, end_id, end_layer, UAV_type, priority;
    char est[10], lat[10];

    while (fscanf(file, "%d %d %d %d %9s %9s %d %d",
                  &start_id, &start_layer, &end_id, &end_layer, est, lat, &UAV_type, &priority) == 8)
    {

        // Belleği yeniden ayarla
        demands = realloc(demands, (count + 1) * sizeof(OD_Demand));
        if (demands == NULL)
        {
            perror("Bellek tahsisi hatası");
            fclose(file);
            return NULL;
        }

        // Veriyi doldur
        demands[count].source_node.node_id = start_id;
        demands[count].source_node.layer = start_layer;
        demands[count].destination_node.node_id = end_id;
        demands[count].destination_node.layer = end_layer;

        // EST ve LAT için belleği ayır
        demands[count].EST = strdup(est);
        // demands[count].ST = strdup(est); // ST is initially the same as EST
        demands[count].LAT = strdup(lat);

        demands[count].UAV_type = UAV_type;
        demands[count].priority = priority;

        (count)++; // OD_Demand sayısını artır
    }

    fclose(file);
    number_of_demand = count;

    demands = realloc(demands, (count + 1) * sizeof(OD_Demand));
    demands[count].source_node.node_id = -1; // Güvenlik önlemi için son elemanı -1 yap

    return demands;
}

double compute_demand_fitness_contribution(String_Rep solution, int od_index)
{
    double f1 = 0, f2 = 0, penalty = 0;

    // f1: path length / max_length
    int od = solution.od_sequence[od_index] - 1; // -1 to match the array index
    int selected = solution.selected_alternative[od_index];
    double length = 0;
    double max_length = find_upperbound(od);

    for (int j = 0; j < NUMBER_OF_EDGES; j++)
    {
        if (alternative_paths_for_ODs[od][selected][j] == 0)
            break;

        int edge = alternative_paths_for_ODs[od][selected][j];
        length += UATFM_route_network[edge].length;
    }

    f1 = length / max_length;
    // tekrar number of paths involveda bölmeye gerek yok zaten bu bölüm herkes için eşit olacağı
    // için etkisi olmaz

    // f2: delay ratio
    if (solution.max_delays[od_index] != 0)
    {
        double delay_ratio = (double)solution.real_delays[od_index] / solution.max_delays[od_index];
        f2 = delay_ratio;
        // tekrar number of paths involveda bölmeye gerek yok zaten bu bölüm herkes için eşit olacağı
        // için etkisi olmaz
    }

    // penalty
    if (solution.real_delays[od_index] > solution.max_delays[od_index])
    {
        penalty = 100000.0; // bu OD gecikmeyi aştıysa tüm çözümün penaltısı onun yüzünden
    }

    double fitness_contrib = (W1 * f1) + (W2 * f2) + penalty;

    return fitness_contrib;
}

// Function to calculate fitness for the entire population
void evaluate_population_fitness(String_Rep *population, double *fitness_values)
{
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        fitness_values[i] = compute_fitness(population[i], i);
    }
}

void reset_edge_matrix()
{
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        for (int l = 0; l < NUMBER_OF_TIME_INTERVALS; l++)
        {
            for (int e = 0; e < NUMBER_OF_EDGES; e++)
            {
                memset(edge_matrix[i][l][e], -1, NUMBER_OF_LINES * sizeof(int));
            }
        }
    }
}

void reset_edge_matrix_for_individual(int individual_index)
{
    for (int l = 0; l < NUMBER_OF_TIME_INTERVALS; l++)
    {
        for (int e = 0; e < NUMBER_OF_EDGES; e++)
        {
            memset(edge_matrix[individual_index][l][e], -1, NUMBER_OF_LINES * sizeof(int));
        }
    }
}

void fill_edge_matrix(String_Rep *population)
{
    int individual_length = 0;

    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        individual_length = get_string_length(population[i]);

        for (int j = 0; j < individual_length; j++)
        {
            int od = population[i].od_sequence[j] - 1; // -1 to match the array index
            int selected = population[i].selected_alternative[j];
            char *EST = demands[od].EST;
            char EST_for_edge_matrix[9];

            subtract_times(EST, "08:00:00", EST_for_edge_matrix); // EST - 08:00:00 == "HH:MM:SS"
            int start_zone = string_time_to_seconds(EST_for_edge_matrix) / 30;

            if (start_zone < 0)
            {
                fprintf(stderr, "[WARN] Negative start_zone corrected from %d to 0\n", start_zone);
                start_zone = 0;
            }
            int initial_start_zone = start_zone;

            double flight_time;
            int selected_edge;

            population[i].real_delays[j] = 0; // real_delay sıfırlanıyor tekrar hesaplama için.

            for (int k = 0; k < NUMBER_OF_EDGES; k++)
            {
                selected_edge = alternative_paths_for_ODs[od][selected][k];

                if (selected_edge == 0) // edge yoksa döngüyü kır
                    break;

                if (k == 0)
                {
                    flight_time = 0;
                    continue;
                }
                else
                {
                    flight_time = UATFM_route_network[selected_edge].length / demands[od].UAV_type; // km / km/h = h
                    flight_time *= 3600;                                                            // saniyeye çevir
                }

                int number_of_zones_to_be_filled = (int)ceil(flight_time / 30.0);

                int number_of_avaiable_zones = 0;
                int starting_point = 0;
                int start_zones[NUMBER_OF_LANES];

                for (int m = 0; m < NUMBER_OF_LANES; m++)
                {
                    start_zones[m] = INT_MAX; // Başlangıç için çok büyük değer
                    starting_point = initial_start_zone;
                    number_of_avaiable_zones = 0;

                    for (int l = initial_start_zone; l < NUMBER_OF_TIME_INTERVALS; l++)
                    {
                        if (number_of_avaiable_zones == number_of_zones_to_be_filled)
                        {
                            int candidate_start = starting_point - number_of_zones_to_be_filled;
                            if (candidate_start < 0)
                                candidate_start = 0;
                            int start_zone_m = candidate_start;
                            int end_zone = start_zone_m + number_of_avaiable_zones;

                            start_zones[m] = start_zone_m;
                            break;
                        }

                        if (edge_matrix[i][l][selected_edge][m] == -1)
                        {
                            number_of_avaiable_zones++;
                        }
                        else
                        {
                            number_of_avaiable_zones = 0;
                        }
                        starting_point++;
                    }
                }

                // En küçük start_zone'u bul
                int min = start_zones[0];
                int min_index = 0;
                for (int r = 1; r < NUMBER_OF_LANES; r++)
                {
                    if (start_zones[r] < min)
                    {
                        min = start_zones[r];
                        min_index = r;
                    }
                }

                start_zone = start_zones[min_index];

                // Eğer hala INT_MAX ise, yani uygun zaman yoksa
                if (start_zone == INT_MAX)
                {
                    fprintf(stderr, "[WARN] No available start zone found for individual %d, od %d, alternative %d\n", i, od, selected);
                    start_zone = initial_start_zone; // ya da başka mantıklı bir default
                }

                int end_zone = start_zone + number_of_zones_to_be_filled;

                // edge_matrix'i doldur
                for (int e = start_zone; e < end_zone && e < NUMBER_OF_TIME_INTERVALS; e++)
                {
                    edge_matrix[i][e][selected_edge][min_index] = od;
                }

                int delay_in_sec = (start_zone - initial_start_zone) * 30;
                population[i].real_delays[j] = delay_in_sec;
            }
        }
    }
}

void fill_edge_matrix_for_individual_with_stop_index(String_Rep *rep, int stop_index, int population_index)
{
    int individual_length = get_string_length(*rep);
    for (int j = 0; j < individual_length; j++)
    {

        if (rep->od_sequence[j] != -1)
        {

            int od = rep->od_sequence[j] - 1; // -1 to match the array index
            int selected = rep->selected_alternative[j];
            // char *ST = demands[od].ST;
            char *EST = demands[od].EST;
            char start_time_for_edge_matrix[9];
            char EST_for_edge_matrix[9];

            subtract_times(EST, "08:00:00", EST_for_edge_matrix); // EST - 08:00:00 == "HH:MM:SS"
            // int EST_zone = string_time_to_seconds(EST_for_edge_matrix) / 30;

            // subtract_times(ST, "08:00:00", start_time_for_edge_matrix); // ST - 08:00:00 == "HH:MM:SS"
            int start_zone = string_time_to_seconds(EST_for_edge_matrix) / 30;
            int initial_start_zone = start_zone;
            double flight_time;
            int selected_edge;

            rep->real_delays[j] = 0; // tekrar hesaplayacağımız için real delayi sıfırlıyoruz.

            for (int k = 0; k < NUMBER_OF_EDGES; k++) // edge sayısı 2048
            {
                selected_edge = alternative_paths_for_ODs[od][selected][k];

                // ODnin seçili alternatifindeki edge'ler bittiğinde -> break
                if (selected_edge == 0)
                {
                    break;
                }

                if (k == 0)
                {
                    flight_time = 0;
                    continue;
                }
                else
                {
                    flight_time = UATFM_route_network[selected_edge].length / demands[od].UAV_type; // km / km/h = h
                    flight_time *= 3600;
                }

                int number_of_zones_to_be_filled = ceil(flight_time / 30);
                int end_zone = start_zone + number_of_zones_to_be_filled;

                int number_of_avaiable_zones;
                int starting_point = 0;
                int start_zones[NUMBER_OF_LANES];

                for (int m = 0; m < NUMBER_OF_LANES; m++)
                {
                    starting_point = initial_start_zone;
                    number_of_avaiable_zones = 0;
                    for (int l = initial_start_zone; l < NUMBER_OF_TIME_INTERVALS; l++)
                    {

                        if (number_of_avaiable_zones == number_of_zones_to_be_filled)
                        {
                            starting_point = starting_point - number_of_zones_to_be_filled;
                            start_zone = starting_point;
                            end_zone = starting_point + number_of_avaiable_zones;
                            start_zones[m] = start_zone;
                            break;
                        }

                        if (edge_matrix[population_index][l][selected_edge][m] == -1)
                        {
                            number_of_avaiable_zones++;
                        }
                        else
                        {
                            number_of_avaiable_zones = 0;
                        }
                        starting_point++;
                    }
                }

                int min = start_zones[0]; // Dizinin ilk elemanını varsayalım ki min olsun
                int min_index = 0;

                for (int r = 1; r < NUMBER_OF_LANES; r++)
                {
                    if (start_zones[r] < min)
                    {
                        min = start_zones[r];
                        min_index = r;
                    }
                }

                start_zone = start_zones[min_index];
                end_zone = start_zone + number_of_avaiable_zones;

                for (int e = start_zone; e < end_zone; e++)
                {
                    edge_matrix[population_index][e][selected_edge][min_index] = od;
                }

                int delay_in_sec = (start_zone - initial_start_zone) * 30;

                rep->real_delays[j] = delay_in_sec;
            }
        }
    }
}

void fill_edge_matrix_for_individual(String_Rep *rep, int individual_index)
{
    int individual_length = get_string_length(*rep);

    for (int j = 0; j < individual_length; j++)
    {

        int od = rep->od_sequence[j] - 1; // -1 to match the array index
        int selected = rep->selected_alternative[j];
        // char *ST = demands[od].ST;
        char *EST = demands[od].EST;
        char start_time_for_edge_matrix[9];
        char EST_for_edge_matrix[9];

        subtract_times(EST, "08:00:00", EST_for_edge_matrix); // EST - 08:00:00 == "HH:MM:SS"
        // int EST_zone = string_time_to_seconds(EST_for_edge_matrix) / 30;

        // subtract_times(ST, "08:00:00", start_time_for_edge_matrix); // ST - 08:00:00 == "HH:MM:SS"
        int start_zone = string_time_to_seconds(EST_for_edge_matrix) / 30;
        int initial_start_zone = start_zone;
        double flight_time;
        int selected_edge;

        rep->real_delays[j] = 0; // tekrar hesaplayacağımız için real delayi sıfırlıyoruz.

        for (int k = 0; k < NUMBER_OF_EDGES; k++) // edge sayısı 2048
        {
            selected_edge = alternative_paths_for_ODs[od][selected][k];

            // ODnin seçili alternatifindeki edge'ler bittiğinde -> break
            if (selected_edge == 0)
            {
                break;
            }

            if (k == 0)
            {
                flight_time = 0;
                continue;
            }
            else
            {
                flight_time = UATFM_route_network[selected_edge].length / demands[od].UAV_type; // km / km/h = h
                flight_time *= 3600;
            }

            int number_of_zones_to_be_filled = ceil(flight_time / 30);
            int end_zone = start_zone + number_of_zones_to_be_filled;

            int number_of_avaiable_zones;
            int starting_point = 0;
            int start_zones[NUMBER_OF_LANES];

            for (int m = 0; m < NUMBER_OF_LANES; m++)
            {
                starting_point = initial_start_zone;
                number_of_avaiable_zones = 0;
                for (int l = initial_start_zone; l < NUMBER_OF_TIME_INTERVALS; l++)
                {

                    if (number_of_avaiable_zones == number_of_zones_to_be_filled)
                    {
                        starting_point = starting_point - number_of_zones_to_be_filled;
                        start_zone = starting_point;
                        end_zone = starting_point + number_of_avaiable_zones;
                        start_zones[m] = start_zone;
                        break;
                    }

                    if (edge_matrix[individual_index][l][selected_edge][m] == -1)
                    {
                        number_of_avaiable_zones++;
                    }
                    else
                    {
                        number_of_avaiable_zones = 0;
                    }
                    starting_point++;
                }
            }

            int min = start_zones[0]; // Dizinin ilk elemanını varsayalım ki min olsun
            int min_index = 0;

            for (int r = 1; r < NUMBER_OF_LANES; r++)
            {
                if (start_zones[r] < min)
                {
                    min = start_zones[r];
                    min_index = r;
                }
            }

            start_zone = start_zones[min_index];
            end_zone = start_zone + number_of_avaiable_zones;

            for (int e = start_zone; e < end_zone; e++)
            {
                edge_matrix[individual_index][e][selected_edge][min_index] = od;
            }

            int delay_in_sec = (start_zone - initial_start_zone) * 30;

            rep->real_delays[j] = delay_in_sec;
        }
    }
}

int find_worst_individual_index(double *fitness_values)
{
    double worst_fitness = DBL_MIN;
    int worst_index = -1;

    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        if (fitness_values[i] > worst_fitness)
        {
            worst_fitness = fitness_values[i];
            worst_index = i;
        }
    }
    return worst_index;
}

void type1_mutation_update_alternatif_path(String_Rep *population, int individual_index)
{
    double randd = rand();
    int length = get_string_length(population[individual_index]);

    if (length == 0)
    {
        return; // Eğer çözüm boşsa devam et
    }

    int rand_od = rand() % length; // Rastgele bir OD seç
    char previous_alt = population[individual_index].selected_alternative[rand_od];

    char new_alt;
    do
    {
        new_alt = rand() % NUMBER_OF_ALTERNATIVES;
    } while (new_alt == previous_alt);

    population[individual_index].selected_alternative[rand_od] = new_alt;

    reset_edge_matrix_for_individual(individual_index);
    fill_edge_matrix_for_individual(&population[individual_index], individual_index);
}

void type3_mutation_add_unplaced_OD(String_Rep *population, int individual_index)
{
    int unplaced_index = -1;
    int individual_length = get_string_length(population[individual_index]);

    // 1. MAX_DELAY'i aşan bir OD bul
    for (int j = 0; j < individual_length; j++)
    {

        if (population[individual_index].real_delays[j] > population[individual_index].max_delays[j])
        {
            // Bu OD'nin indexini al
            unplaced_index = j;
            break;
        }
    }

    if (unplaced_index == -1)
        return;

    int unplaced_demand = population[individual_index].od_sequence[unplaced_index];
    int unplaced_demand_priority = demands[unplaced_demand].priority;

    // 2. Aynı priority'e sahip diğer OD indexlerini topla (unplaced dışında)
    int matching_demands[individual_length];
    int match_count = 0;

    for (int k = 0; k < individual_length; k++)
    {
        if (k == unplaced_index)
            continue;

        int od2 = population[individual_index].od_sequence[k];
        if (demands[od2].priority == unplaced_demand_priority)
        {
            matching_demands[match_count++] = k;
        }
    }

    if (match_count == 0)
        return;

    // 3. Random olarak bir tanesini seç
    int rand_idx = rand() % match_count;
    int swap_index = matching_demands[rand_idx];

    // 4. Swap işlemi
    int temp = population[individual_index].od_sequence[unplaced_index];
    population[individual_index].od_sequence[unplaced_index] = population[individual_index].od_sequence[swap_index];
    population[individual_index].od_sequence[swap_index] = temp;

    reset_edge_matrix_for_individual(individual_index);
    fill_edge_matrix_for_individual(&population[individual_index], individual_index);
}

void type4_mutation_change_OD_order(String_Rep *population, int individual_index)
{
    double randd = rand();
    int length = get_string_length(population[individual_index]);

    if (length == 0)
    {
        return; // Eğer çözüm boşsa devam et
    }

    int index1 = rand() % length; // Rastgele bir OD seç

    int od1 = population[individual_index].od_sequence[index1];
    int priority1 = demands[od1].priority;

    int index2;
    int attempts = 0;

    do
    {
        int index2 = rand() % length; // Aynı önceliğe sahip başka bir OD seç
        int od2 = population[individual_index].od_sequence[index2];

        if (index1 != index2 && demands[od2].priority == priority1)
        {
            // Sırayı değiştir
            int temp = population[individual_index].od_sequence[index1];
            population[individual_index].od_sequence[index1] = population[individual_index].od_sequence[index2];
            population[individual_index].od_sequence[index2] = temp;

            break;
        }

        attempts++;

    } while (attempts < length / 3); // Sonsuz döngüden kaçınmak için deneme sınırı

    reset_edge_matrix_for_individual(individual_index);
    fill_edge_matrix_for_individual(&population[individual_index], individual_index);
}

void type5_mutation_update_alternatif_path_based_on_fitness(String_Rep *population, int individual_index)
{
    // 1. En kötü katkıyı yapan OD’yi bul
    double worst_contrib = -1.0;
    int worst_od_index = -1;
    int individual_length = get_string_length(population[individual_index]);

    for (int j = 0; j < individual_length; j++)
    {
        double contrib = compute_demand_fitness_contribution(population[individual_index], j);
        if (contrib > worst_contrib)
        {
            worst_contrib = contrib;
            worst_od_index = j;
        }
    }

    if (worst_od_index == -1)
        return;

    int current_alt = population[individual_index].selected_alternative[worst_od_index];
    int best_alt = current_alt;
    double best_contrib = worst_contrib;

    // 2. Alternatif path’leri dene
    for (int alt = 0; alt < NUMBER_OF_ALTERNATIVES; alt++)
    {
        if (alt == current_alt)
            continue;

        // Geçici çözüm oluştur
        String_Rep temp = population[individual_index];
        temp.selected_alternative[worst_od_index] = alt;

        // Yeni katkıyı hesapla
        double new_contrib = compute_demand_fitness_contribution(temp, worst_od_index);

        if (new_contrib < best_contrib)
        {
            best_contrib = new_contrib;
            best_alt = alt;
        }
    }

    // 3. En iyi alternatifi uygula
    if (best_alt != current_alt)
    {
        population[individual_index].selected_alternative[worst_od_index] = best_alt;
        reset_edge_matrix_for_individual(individual_index);
        fill_edge_matrix_for_individual(&population[individual_index], individual_index);
    }
}
// GAA
void generate_offspring(String_Rep *population, double *fitness_values, int decision, int parent1_index, int parent2_index, int worst_index)
{
    String_Rep child1;
    String_Rep child2;

    for (int i = 0; i < number_of_demand + 1; i++)
    {
        child1.od_sequence[i] = 0;
        child2.od_sequence[i] = 0;
    }
    String_Rep parent1 = population[parent1_index];
    String_Rep parent2 = population[parent2_index];

    if (decision == 0)
    {
        order_crossover(parent1, parent2, &child1, &child2);

        memcpy(&population[worst_index], &child1, sizeof(String_Rep));
        reset_edge_matrix_for_individual(worst_index);
        fill_edge_matrix_for_individual(&population[worst_index], worst_index);

        int child1_fitness = compute_fitness(population[worst_index], worst_index);

        memcpy(&population[worst_index], &child2, sizeof(String_Rep));
        reset_edge_matrix_for_individual(worst_index);
        fill_edge_matrix_for_individual(&population[worst_index], worst_index);

        int child2_fitness = compute_fitness(population[worst_index], worst_index);

        if (child1_fitness > child2_fitness)
        {
            memcpy(&child1, &child2, sizeof(String_Rep));
        }
    }

    else
    {
        if (find_similarity_ratio(parent1, parent2) > 80.0f)
        {
            if (decision == 1)
            {
                similarity_crossover_minimize_conflict(parent1, parent2, &child1);
            }
            else
            {
                similarity_crossover_fitness_related(parent1, parent2, &child1);
            }
        }
        else
        {
            // Kopyalama
            memcpy(edge_matrix_for_individual1, edge_matrix[parent1_index], sizeof(edge_matrix_for_individual1));
            memcpy(edge_matrix_for_individual2, edge_matrix[parent2_index], sizeof(edge_matrix_for_individual2));

            low_similarity_crossover(parent1, parent2, &child1, &child2, population, parent1_index, parent2_index);
            memcpy(edge_matrix[parent1_index], edge_matrix_for_individual1, sizeof(edge_matrix_for_individual1));
            memcpy(edge_matrix[parent2_index], edge_matrix_for_individual2, sizeof(edge_matrix_for_individual2));

            memcpy(&population[worst_index], &child1, sizeof(String_Rep));
            reset_edge_matrix_for_individual(worst_index);
            fill_edge_matrix_for_individual(&population[worst_index], worst_index);

            int child1_fitness = compute_fitness(population[worst_index], worst_index);

            memcpy(&population[worst_index], &child2, sizeof(String_Rep));
            reset_edge_matrix_for_individual(worst_index);
            fill_edge_matrix_for_individual(&population[worst_index], worst_index);

            int child2_fitness = compute_fitness(population[worst_index], worst_index);

            if (child1_fitness > child2_fitness)
            {
                memcpy(&child1, &child2, sizeof(String_Rep));
            }
        }
    }

    memcpy(&population[worst_index], &child1, sizeof(String_Rep));

    // Reinitialize the edge matrix for the new population
    reset_edge_matrix_for_individual(worst_index);
    fill_edge_matrix_for_individual(&population[worst_index], worst_index);
}

void copy_parent(String_Rep *source, String_Rep *destination)
{
    // od_sequence, selected_alternative, max_delays ve real_delays dizilerini kopyala
    int length = get_string_length(*source);
    for (int i = 0; i < length; i++)
    {
        destination->od_sequence[i] = source->od_sequence[i];
        destination->selected_alternative[i] = source->selected_alternative[i];
        destination->max_delays[i] = source->max_delays[i];
        destination->real_delays[i] = source->real_delays[i];
    }

    if (length != NUMBER_OF_ODS * MULTIPLICITY_STRING_REP)
    {
        destination->od_sequence[length] = 0;
    }
}
void try_unplaced_ODs_to_put_solution(int unplaced_ODs[], String_Rep *child, int replication_index, int treshold, int length, int population_index)
{
    for (int i = 0; i < number_of_demand; i++)
    {

        if (unplaced_ODs[i] == 0)
        {
            child->od_sequence[replication_index] = -1;
            child->selected_alternative[replication_index] = 0;
            child->max_delays[replication_index] = 0;
            child->real_delays[replication_index] = 0;

            break;
        }

        else
        {

            if (unplaced_ODs[i] != -1)
            {
                short temp_od = unplaced_ODs[i];
                child->od_sequence[replication_index] = temp_od;

                char temp_alt = 0;
                child->selected_alternative[replication_index] = temp_alt;

                reset_edge_matrix_for_individual(population_index);
                fill_edge_matrix_for_individual_with_stop_index(child, replication_index, population_index);

                child->max_delays[replication_index] = find_max_delay(*child, replication_index);
                if (child->real_delays[replication_index] - child->max_delays[replication_index] < treshold)
                {
                    unplaced_ODs[i] = -1;
                    reset_edge_matrix_for_individual(population_index);
                    break;
                }

                reset_edge_matrix_for_individual(population_index);
            }
        }
    }
}

void calculate_unplaced_ODs(String_Rep child, int not_found_ods[], int length)
{
    int od = 1;
    int count = 0;

    for (int m = 0; m < number_of_demand; m++)
    {
        int od_is_found = 0;

        // TODO: number od demands yerine uzunluk fonksiyonu çağırılacak
        for (int i = 0; i < length; i++)
        {
            if (child.od_sequence[i] == od)
            {
                od_is_found = 1;
                break;
            }
        }

        if (od_is_found == 0)
        {
            not_found_ods[count] = od;
            count++;
        }

        od++;
    }
}

void calculate_repeated_ods_indexes(String_Rep child, int repeated_od_indexes[], int length)
{
    String_Rep childCopy;
    memset(&childCopy, 0, sizeof(String_Rep));
    copy_parent(&child, &childCopy);

    int index = 0;

    // buraya da uzunluk fonksiyonu çağırılıcak
    for (int i = 0; i < length - 2; i++)
    {
        if (childCopy.od_sequence[i] == 0)
        {
            break;
        }

        for (int j = i + 1; j < length - 1; j++)
        {

            if (childCopy.od_sequence[i] == -1)
            {
                break;
            }

            if (childCopy.od_sequence[i] == childCopy.od_sequence[j])
            {
                repeated_od_indexes[index++] = j;
                childCopy.od_sequence[j] = -1;
            }
        }
    }
}

void low_similarity_crossover(String_Rep parent1, String_Rep parent2, String_Rep *child1, String_Rep *child2, String_Rep *population, int parent1_index, int parent2_index)
{
    int parent1_length = get_string_length(parent1);
    int parent2_length = get_string_length(parent2);

    int max_length = MAX(parent1_length, parent2_length);

    int maxInterval = ceil(sqrt(number_of_demand));

    int isRealDelaysGreaterThanMaxDealys[max_length][2];

    for (int i = 0; i < max_length; i++)
    {

        isRealDelaysGreaterThanMaxDealys[i][0] = 0;
        isRealDelaysGreaterThanMaxDealys[i][1] = 0;
    }

    String_Rep parent1Copy, parent2Copy;

    int numberOfIterationsForParent1 = floor(get_string_length(parent1) - maxInterval) + 1;
    int numberOfIterationsForParent2 = floor(get_string_length(parent2) - maxInterval) + 1;

    int start = 0;
    int end = start + maxInterval;
    int maxCountForParent1 = 0;
    int maxCountForParent2 = 0;
    int bestStartForParent1 = 0;
    int bestEndForParent1 = start + maxInterval;
    int bestStartForParent2 = 0;
    int bestEndForParent2 = start + maxInterval;
    int countForParent1 = 0;
    int countForParent2 = 0;
    int treshold = 30;

    for (int i = 0; i < max_length; i++)
    {

        if (parent1.od_sequence[i] != 0)
        {
            if (parent1.real_delays[i] > parent1.max_delays[i])
            {
                isRealDelaysGreaterThanMaxDealys[i][0] = 1;
            }
            else
            {
                isRealDelaysGreaterThanMaxDealys[i][0] = 0;
            }
        }
        else
        {
            break;
        }
    }

    for (int i = 0; i < max_length; i++)
    {
        if (parent2.od_sequence[i] != 0)
        {
            if (parent2.real_delays[i] > parent2.max_delays[i])
            {
                isRealDelaysGreaterThanMaxDealys[i][1] = 1;
            }
            else
            {
                isRealDelaysGreaterThanMaxDealys[i][1] = 0;
            }
        }
        else
        {
            break;
        }
    }

    for (int i = 0; i < numberOfIterationsForParent1; i++)
    {
        countForParent1 = 0;

        for (int j = start; j < end; j++)
        {

            if (isRealDelaysGreaterThanMaxDealys[j][0] == 1)
            {
                countForParent1++;
            }
        }

        if (countForParent1 > maxCountForParent1)
        {
            maxCountForParent1 = countForParent1;
            bestStartForParent1 = start;
            bestEndForParent1 = end;
        }

        start = start + 1;
        end = start + maxInterval;
    }

    start = 0;
    end = start + maxInterval;

    for (int i = 0; i < numberOfIterationsForParent2; i++)
    {
        countForParent2 = 0;

        for (int j = start; j < end; j++)
        {

            if (isRealDelaysGreaterThanMaxDealys[j][1] == 1)
            {
                countForParent2++;
            }
        }

        if (countForParent2 > maxCountForParent2)
        {
            maxCountForParent2 = countForParent2;
            bestStartForParent2 = start;
            bestEndForParent2 = end;
        }

        start = start + 1;
        end = start + maxInterval;
    }

    copy_parent(&parent1, &parent1Copy);
    copy_parent(&parent2, &parent2Copy);

    for (int a = bestStartForParent1; a < bestEndForParent1; a++)
    {
        parent1Copy.max_delays[a] = 0;
        parent1Copy.od_sequence[a] = -1;
        parent1Copy.selected_alternative[a] = 0;
        parent1Copy.real_delays[a] = 0;
    }

    for (int a = bestStartForParent2; a < bestEndForParent2; a++)
    {
        parent2Copy.max_delays[a] = 0;
        parent2Copy.od_sequence[a] = -1;
        parent2Copy.selected_alternative[a] = 0;
        parent2Copy.real_delays[a] = 0;
    }

    int m = bestStartForParent2;

    for (int j = bestStartForParent1; j < bestEndForParent1; j++)
    {

        if (parent1.real_delays[j] - parent1.max_delays[j] < treshold && parent2.real_delays[m] - parent2.max_delays[m] > treshold)
        {

            short temp_od = parent1.od_sequence[j];
            parent1Copy.od_sequence[j] = temp_od;
            parent2Copy.od_sequence[m] = temp_od;

            char temp_alt = parent1.selected_alternative[j];
            parent1Copy.selected_alternative[j] = temp_alt;
            parent2Copy.selected_alternative[m] = temp_alt;

            int temp_maxdelay = parent1.max_delays[j];
            parent1Copy.max_delays[j] = temp_maxdelay;
            parent2Copy.max_delays[m] = temp_maxdelay;
        }

        if (parent1.real_delays[j] - parent1.max_delays[j] > treshold && parent2.real_delays[m] - parent2.max_delays[m] < treshold)
        {
            short temp_od = parent2.od_sequence[m];
            parent1Copy.od_sequence[j] = temp_od;
            parent2Copy.od_sequence[m] = temp_od;

            char temp_alt = parent2.selected_alternative[m];
            parent1Copy.selected_alternative[j] = temp_alt;
            parent2Copy.selected_alternative[m] = temp_alt;

            int temp_maxdelay = parent2.max_delays[m];
            parent1Copy.max_delays[j] = temp_maxdelay;
            parent2Copy.max_delays[m] = temp_maxdelay;
        }

        if (parent1.real_delays[j] - parent1.max_delays[j] < treshold && parent2.real_delays[m] - parent2.max_delays[m] < treshold)
        {
            int delay1 = parent1.real_delays[j] - parent1.max_delays[j];
            int delay2 = parent2.real_delays[m] - parent2.max_delays[m];

            if (delay1 <= delay2)
            {
                short temp_od = parent1.od_sequence[j];
                parent1Copy.od_sequence[j] = temp_od;
                parent2Copy.od_sequence[m] = temp_od;

                char temp_alt = parent1.selected_alternative[j];
                parent1Copy.selected_alternative[j] = temp_alt;
                parent2Copy.selected_alternative[m] = temp_alt;

                int temp_maxdelay = parent1.max_delays[j];
                parent1Copy.max_delays[j] = temp_maxdelay;
                parent2Copy.max_delays[m] = temp_maxdelay;
            }

            if (delay2 < delay1)
            {
                short temp_od = parent2.od_sequence[m];
                parent1Copy.od_sequence[j] = temp_od;
                parent2Copy.od_sequence[m] = temp_od;

                char temp_alt = parent2.selected_alternative[m];
                parent1Copy.selected_alternative[j] = temp_alt;
                parent2Copy.selected_alternative[m] = temp_alt;

                int temp_maxdelay = parent2.max_delays[m];
                parent1Copy.max_delays[j] = temp_maxdelay;
                parent2Copy.max_delays[m] = temp_maxdelay;
            }
        }

        m++;
    }

    m = bestStartForParent2;

    for (int j = bestStartForParent1; j < bestEndForParent1; j++)
    {

        if (parent1.real_delays[j] - parent1.max_delays[j] > treshold && parent2.real_delays[m] - parent2.max_delays[m] > treshold)
        {

            int unplacedODsForChild1[number_of_demand];
            int unplacedODsForChild2[number_of_demand];

            for (int i = 0; i < number_of_demand; i++)
            {
                unplacedODsForChild1[i] = 0;
                unplacedODsForChild2[i] = 0;
            }

            calculate_unplaced_ODs(parent1Copy, unplacedODsForChild1, parent1_length);
            calculate_unplaced_ODs(parent2Copy, unplacedODsForChild2, parent2_length);

            try_unplaced_ODs_to_put_solution(unplacedODsForChild1, &parent1Copy, j, treshold, parent1_length, parent1_index);
            try_unplaced_ODs_to_put_solution(unplacedODsForChild2, &parent2Copy, m, treshold, parent2_length, parent2_index);
        }
        m++;
    }

    cancel_invalid_demands(&parent1Copy, parent1_length);
    cancel_invalid_demands(&parent2Copy, parent2_length);

    int findReapetedOdIndexesForChild1[max_length];
    int findReapetedOdIndexesForChild2[max_length];
    for (int i = 0; i < max_length; i++)
    {
        findReapetedOdIndexesForChild1[i] = 0;
        findReapetedOdIndexesForChild2[i] = 0;
    }
    calculate_repeated_ods_indexes(parent1Copy, findReapetedOdIndexesForChild1, parent1_length);
    calculate_repeated_ods_indexes(parent2Copy, findReapetedOdIndexesForChild2, parent2_length);

    // Print the lengths of findReapetedOdIndexesForChild1 and findReapetedOdIndexesForChild2
    int len1 = 0, len2 = 0;
    for (int i = 0; i < max_length; i++)
    {
        if (findReapetedOdIndexesForChild1[i] != 0)
            len1++;
        if (findReapetedOdIndexesForChild2[i] != 0)
            len2++;
    }
    int unplacedOdsOfChild1[number_of_demand];
    int unplacedOdsOfChild2[number_of_demand];

    for (int i = 0; i < number_of_demand; i++)
    {
        unplacedOdsOfChild1[i] = 0;
        unplacedOdsOfChild2[i] = 0;
    }

    calculate_unplaced_ODs(parent1Copy, unplacedOdsOfChild1, parent1_length);
    calculate_unplaced_ODs(parent2Copy, unplacedOdsOfChild2, parent2_length);

    // Print the lengths of unplacedOdsOfChild1 and unplacedOdsOfChild2
    int unplaced1_len = 0, unplaced2_len = 0;
    for (int i = 0; i < number_of_demand; i++)
    {
        if (unplacedOdsOfChild1[i] != 0)
            unplaced1_len++;
        if (unplacedOdsOfChild2[i] != 0)
            unplaced2_len++;
    }

    for (int a = 0; a < max_length; a++)
    {

        int index = findReapetedOdIndexesForChild1[a];

        if (index == 0)
        {
            break;
        }

        try_unplaced_ODs_to_put_solution(unplacedOdsOfChild1, &parent1Copy, index, treshold, parent1_length, parent1_index);
        findReapetedOdIndexesForChild1[a] = 0;
    }

    for (int a = 0; a < max_length; a++)
    {

        int index = findReapetedOdIndexesForChild2[a];

        if (index == 0)
        {
            break;
        }

        try_unplaced_ODs_to_put_solution(unplacedOdsOfChild2, &parent2Copy, index, treshold, parent2_length, parent2_index);
        findReapetedOdIndexesForChild2[a] = 0;
    }

    cancel_invalid_demands(&parent1Copy, parent1_length);
    cancel_invalid_demands(&parent2Copy, parent2_length);

    reset_edge_matrix_for_individual(parent1_index);
    reset_edge_matrix_for_individual(parent2_index);

    int parent1_cpy_length = get_string_length(parent1Copy);
    int parent2_cpy_length = get_string_length(parent2Copy);

    for (int i = 0; i < parent1_cpy_length; i++)
    {
        child1->od_sequence[i] = parent1Copy.od_sequence[i];
        child1->selected_alternative[i] = parent1Copy.selected_alternative[i];
        child1->max_delays[i] = parent1Copy.max_delays[i];
        child1->real_delays[i] = parent1Copy.real_delays[i];
    }
    for (int i = 0; i < parent2_cpy_length; i++)
    {
        child2->od_sequence[i] = parent2Copy.od_sequence[i];
        child2->selected_alternative[i] = parent2Copy.selected_alternative[i];
        child2->max_delays[i] = parent2Copy.max_delays[i];
        child2->real_delays[i] = parent2Copy.real_delays[i];
    }

    if (parent1_cpy_length != NUMBER_OF_ODS * MULTIPLICITY_STRING_REP)
    {
        child1->od_sequence[parent1_cpy_length] = 0;
    }

    if (parent2_cpy_length != NUMBER_OF_ODS * MULTIPLICITY_STRING_REP)
    {
        child2->od_sequence[parent2_cpy_length] = 0;
    }
}

/*---------------------------------------------------------------------*/

void print_best(String_Rep best_solution)
{
    printf("OD sequence: ");

    for (int i = 0; i < get_string_length(best_solution); i++)
    {
        printf("%d ", best_solution.od_sequence[i]);
    }

    printf("\n Selected Alternative: ");

    for (int i = 0; i < get_string_length(best_solution); i++)
    {
        printf("%d ", best_solution.selected_alternative[i]);
    }

    printf("\nMax delays: ");
    for (int i = 0; i < get_string_length(best_solution); i++)
    {
        printf("%d ", best_solution.max_delays[i]);
    }

    printf("\n\nReal delays: ");

    for (int i = 0; i < get_string_length(best_solution); i++)
    {
        printf("%d ", best_solution.real_delays[i]);
    }
}

double average_fitness(double *fitness_values)
{
    double sum = 0;
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        sum += fitness_values[i];
    }
    return sum / POPULATION_SIZE;
}


int init_test(const char* edge_file_path, const char* demand_file_path, const char* path_file_path) {
    
    
    if (demands != NULL) {
        for (int i = 0; i < number_of_demand; i++) {
            if (demands[i].EST != NULL) {
                memset(demands[i].EST, 0, MAX_TIME_STRING_LENGTH); // sabit uzunlukla
            }
            if (demands[i].LAT != NULL) {
                memset(demands[i].LAT, 0, MAX_TIME_STRING_LENGTH);
            }
        }
        memset(demands, 0, number_of_demand * sizeof(OD_Demand));
    }

    // 3D diziler - içeriklerini de sıfırla
    for (int i = 0; i < NUMBER_OF_ODS; i++) {
    for (int j = 0; j <= NUMBER_OF_ALTERNATIVES; j++) {
        if (all_alternative_paths_for_ODs[i][j]) {
            memset(all_alternative_paths_for_ODs[i][j], 0, MAX_NUMBER_OF_EDGES_IN_ROUTE * sizeof(unsigned short));
        }
    }
}

    if (alternative_paths_for_ODs != NULL) {
        for (int i = 0; i < NUMBER_OF_ODS * MULTIPLICITY_STRING_REP; i++) {
            for (int j = 0; j < NUMBER_OF_ALTERNATIVES; j++) {
                if (alternative_paths_for_ODs[i][j] != NULL) {
                    memset(alternative_paths_for_ODs[i][j], 0, MAX_NUMBER_OF_EDGES_IN_ROUTE * sizeof(unsigned short));
                }
            }
        }
    }

    // population ve backup sıfırlanabilir çünkü sabit diziler içeriyor
    if (population != NULL) {
        memset(population, 0, POPULATION_SIZE * sizeof(String_Rep));
    }
    if (population_backup != NULL) {
        memset(population_backup, 0, POPULATION_SIZE * sizeof(String_Rep));
    }


    // --- 1. EDGE DOSYASINI OKU ---
    FILE *input_file = fopen(edge_file_path, "r");
    if (input_file == NULL) {
        perror("Error opening edge input file.");
        return 1;
    }

    int edge_id, start_node, end_node, start_layer, end_layer, edge_type;
    float length, total_cost, flight_cost, noise_cost;
    int index = 1;

    while (fscanf(input_file, "%d %d %d %d %d %d %f %f %f %f", &edge_id,
                  &start_node, &end_node, &start_layer, &end_layer, &edge_type,
                  &length, &flight_cost, &noise_cost, &total_cost) == 10)
    {
        UATFM_route_network[index].edge_id = edge_id;
        UATFM_route_network[index].start_node = start_node;
        UATFM_route_network[index].end_node = end_node;
        UATFM_route_network[index].start_layer = start_layer;
        UATFM_route_network[index].end_layer = end_layer;
        UATFM_route_network[index].edge_type = edge_type;
        UATFM_route_network[index].length = length / 1000;
        UATFM_route_network[index].total_cost = total_cost;
        UATFM_route_network[index].flight_cost = flight_cost;
        UATFM_route_network[index].noise_cost = noise_cost;
        index++;
    }

    fclose(input_file);
    
    


    // --- 2. DEMAND DOSYASINI OKU ---
    demands = read_od_demand(demand_file_path);
    if (number_of_demand >= NUMBER_OF_ODS * MULTIPLICITY_STRING_REP) {
        fprintf(stderr, "Warning: Number of demands exceeds the maximum limit.\n");
        return EXIT_FAILURE;
    }

    // --- 3. PATH DOSYASINI OKU ---
    FILE *file = fopen(path_file_path, "r");
    if (file == NULL) {
        perror("Error opening path file");
        return EXIT_FAILURE;
    }

    char line[MAX_LINE_LENGTH];
    int od_index = 0, alt_index = 0, blank_line_count = 0;

    while (fgets(line, sizeof(line), file)) {
        if (is_blank_line(line)) {
            blank_line_count++;
            if (blank_line_count == 1) {
                od_index++;
                alt_index = 0;
                blank_line_count = 0;

                if (od_index > NUMBER_OF_ODS) {
                    printf("od_index: %d\n", od_index);
                    fprintf(stderr, "Warning: OD index exceeded limit.\n");
                    break;
                }
            }
            continue;
        }

        blank_line_count = 0;

        if (alt_index >= NUMBER_OF_ALTERNATIVES + 1) {
            fprintf(stderr, "Warning: Too many alternatives for OD Pair %d\n", od_index);
            continue;
        }

        char *token = strtok(line, " ");
        int edge_count = 0;

        while (token != NULL) {
            if (edge_count >= MAX_NUMBER_OF_EDGES_IN_ROUTE) {
                fprintf(stderr, "Warning: Too many edges in alternative %d of OD Pair %d\n", alt_index, od_index);
                break;
            }

            all_alternative_paths_for_ODs[od_index][alt_index][edge_count] = (unsigned short)atoi(token);
            edge_count++;
            token = strtok(NULL, " ");
        }

        alt_index++;
    }

    fclose(file);

    // --- 4. ALTERNATİFLERİ DEMAND'LERE ATA ---
    for (int i = 0; i < number_of_demand; i++) {
        for (int j = 0; j < NUMBER_OF_ODS; j++) {
            if (demands[i].source_node.node_id == all_alternative_paths_for_ODs[j][0][0] &&
                demands[i].source_node.layer == all_alternative_paths_for_ODs[j][0][1] &&
                demands[i].destination_node.node_id == all_alternative_paths_for_ODs[j][0][2] &&
                demands[i].destination_node.layer == all_alternative_paths_for_ODs[j][0][3]) {
                
                for (int k = 0; k < NUMBER_OF_ALTERNATIVES; k++) {
                    for (int l = 0; 1; l++) {
                        if (all_alternative_paths_for_ODs[j][k + 1][l] == 0) {
                            break;
                        }
                        alternative_paths_for_ODs[i][k][l] = all_alternative_paths_for_ODs[j][k + 1][l];
                    }
                }
            }
        }
    }

    // --- 5. POPÜLASYONU BAŞLAT ---
    for (int i = 0; i < POPULATION_SIZE; i++) {
        if (!population_backup[i].od_sequence ||
            !population_backup[i].selected_alternative ||
            !population_backup[i].max_delays ||
            !population_backup[i].real_delays) {
            perror("Memory allocation failed!");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < number_of_demand; j++) {
            population_backup[i].od_sequence[j] = j + 1;
        }

        shuffle_priority_aware(population_backup[i].od_sequence, number_of_demand);

        for (int j = 0; j < number_of_demand; j++) {
            population_backup[i].selected_alternative[j] = rand() % NUMBER_OF_ALTERNATIVES;
            population_backup[i].max_delays[j] = find_max_delay(population_backup[i], j);
            population_backup[i].real_delays[j] = 0;
        }
    }

    return 0;
}


int test(int crossover_type, int mutation_type)
{
    // şimdilik başlangıç saat 8:00 bitiş 18:00 olarak düşünelim
    memcpy(population, population_backup, sizeof(String_Rep) * POPULATION_SIZE);

    double fitness_values[POPULATION_SIZE];

    // Initialize the edge matrix
    reset_edge_matrix();
    fill_edge_matrix(population);

    evaluate_population_fitness(population, fitness_values);
    // find best solution

    int best_index = 0;
    double best_fitness = fitness_values[0];
    for (int i = 0; i < POPULATION_SIZE; i++)
    {

        if (fitness_values[i] < best_fitness)
        {
            best_fitness = fitness_values[i];
            best_index = i;
        }
    }

    int worst_index_foreval = 0;
    double worst_fitness = fitness_values[0];
    for (int i = 0; i < POPULATION_SIZE; i++)
    {

        if (fitness_values[i] > worst_fitness)
        {
            worst_fitness = fitness_values[i];
            worst_index_foreval = i;
        }
    }
    double population_fittness_average = average_fitness(fitness_values);

    String_Rep best_solution = population[best_index];
    String_Rep worst_solution = population[worst_index_foreval];

    //***** Best Solution Fittness *****

    double *fitness_of_best_solution = compute_fitness_for_best_solution(best_solution, best_index);
    double *fitness_of_worst_solution = compute_fitness_for_best_solution(worst_solution, worst_index_foreval);

    printf("\nSSGA Iteration No: %d F1: %f, F2: %f, F3: %f, Fitness Sum: %f, F1_num: %f, F1_denom: %f, F2_num: %f, F3_num: %f, F3_denom: %f, Penalty Count: %f, Population Fitness Average %f, Worst Population: F1: %f, F2: %f, F3: %f, Fitness Sum: %f, F1_num: %f, F1_denom: %f, F2_num: %f, F3_num: %f, F3_denom: %f, Penalty Count: %f", 0, fitness_of_best_solution[0], fitness_of_best_solution[1], fitness_of_best_solution[2], fitness_of_best_solution[3], fitness_of_best_solution[4], fitness_of_best_solution[5], fitness_of_best_solution[6], fitness_of_best_solution[7], fitness_of_best_solution[8], fitness_of_best_solution[9], population_fittness_average, fitness_of_worst_solution[0], fitness_of_worst_solution[1], fitness_of_worst_solution[2], fitness_of_worst_solution[3], fitness_of_worst_solution[4], fitness_of_worst_solution[5], fitness_of_worst_solution[6], fitness_of_worst_solution[7], fitness_of_worst_solution[8], fitness_of_worst_solution[9]);
    // LOOP FOR GA
    for (int GA_index = 0; GA_index < NUMBER_OF_GENERATIONS; GA_index++)
    {
        // Select parents using tournament selection
        int parent1_index = tournament_selection(population, fitness_values);
        int parent2_index = tournament_selection(population, fitness_values);

        while (parent1_index == parent2_index)
        {
            parent2_index = tournament_selection(population, fitness_values);
        }

        int worst_index = find_worst_individual_index(fitness_values);

        //  generate new population using crossover
        generate_offspring(population, fitness_values, crossover_type, parent1_index, parent2_index, worst_index);
        // 0 ise order crossover, 1 - 2 ise problem specific crossover
        // similarity yüksekse --> 1 ise similarity crossover minimize conflict, 2 ise similarity crossover fitness related

        // mutation
        if ((double)rand() / RAND_MAX < MUTATION_RATE) // MUTATION_RATE ihtimalle mutasyon uygula
        {
            if (mutation_type == 0)
            {
                type1_mutation_update_alternatif_path(population, worst_index);
            }
            else if (mutation_type == 1)
            {
                type3_mutation_add_unplaced_OD(population, worst_index);
            }
            else if (mutation_type == 2)
            {
                type4_mutation_change_OD_order(population, worst_index);
            }
            else if (mutation_type == 3)
            {
                type5_mutation_update_alternatif_path_based_on_fitness(population, worst_index);
            }
        }

        fitness_values[worst_index] = compute_fitness(population[worst_index], worst_index);

        // find best solution

        int best_index = 0;
        double best_fitness = fitness_values[0];
        for (int i = 0; i < POPULATION_SIZE; i++)
        {

            if (fitness_values[i] < best_fitness)
            {
                best_fitness = fitness_values[i];
                best_index = i;
            }
        }

        int worst_index_foreval = 0;
        double worst_fitness = fitness_values[0];
        for (int i = 0; i < POPULATION_SIZE; i++)
        {

            if (fitness_values[i] > worst_fitness)
            {
                worst_fitness = fitness_values[i];
                worst_index_foreval = i;
            }
        }
        double population_fittness_average = average_fitness(fitness_values);

        String_Rep best_solution = population[best_index];
        String_Rep worst_solution = population[worst_index_foreval];

        //***** Best Solution Fittness *****

        double *fitness_of_best_solution = compute_fitness_for_best_solution(best_solution, best_index);
        double *fitness_of_worst_solution = compute_fitness_for_best_solution(worst_solution, worst_index_foreval);

        printf("\nSSGA Iteration No: %d F1: %f, F2: %f, F3: %f, Fitness Sum: %f, F1_num: %f, F1_denom: %f, F2_num: %f, F3_num: %f, F3_denom: %f, Penalty Count: %f, Population Fitness Average %f, Worst Population: F1: %f, F2: %f, F3: %f, Fitness Sum: %f, F1_num: %f, F1_denom: %f, F2_num: %f, F3_num: %f, F3_denom: %f, Penalty Count: %f", GA_index + 1, fitness_of_best_solution[0], fitness_of_best_solution[1], fitness_of_best_solution[2], fitness_of_best_solution[3], fitness_of_best_solution[4], fitness_of_best_solution[5], fitness_of_best_solution[6], fitness_of_best_solution[7], fitness_of_best_solution[8], fitness_of_best_solution[9], population_fittness_average, fitness_of_worst_solution[0], fitness_of_worst_solution[1], fitness_of_worst_solution[2], fitness_of_worst_solution[3], fitness_of_worst_solution[4], fitness_of_worst_solution[5], fitness_of_worst_solution[6], fitness_of_worst_solution[7], fitness_of_worst_solution[8], fitness_of_worst_solution[9]);
    }
}

#ifdef _WIN32
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir) // mode parametresi yok
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

int main()
{
    srand(time(NULL));
    allocate_all();
    time_t rawtime;
    struct tm *timeinfo;
    char folder_name[64];

    // Şu anki zamanı al
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    // outputs_HH_MM_SS formatında klasör ismi oluştur
    strftime(folder_name, sizeof(folder_name), "outputs_%H_%M_%S", timeinfo);

    // mkdir ile klasörü oluştur
    if (mkdir(folder_name, 0777) == 0)
    {
        printf("Klasör oluşturuldu: %s\n", folder_name);
    }
    else if (errno == EEXIST)
    {
        printf("Klasör zaten var: %s\n", folder_name);
    }
    else
    {
        perror("mkdir");
    }


    const char* small_demand_files[] = {
        "demands_small_scale/demands_4_6.txt",
        "demands_small_scale/demands_4_12.txt",
        "demands_small_scale/demands_8_6.txt",
        "demands_small_scale/demands_8_12.txt"
    };

    const char* small_sim_versions[] = {
        "shortest_path_management/formatted_paths/formatted_output_small_sim1.txt",
        "shortest_path_management/formatted_paths/formatted_output_small_sim2.txt"
    };

    const char* small_edge_file = "9138_1243_ist_small_scale.txt";  // sabit veya kullanılmayacaksa NULL olabilir

    for (int i = 0; i < sizeof(small_demand_files) / sizeof(small_demand_files[0]); i++) {
        for (int j = 0; j < sizeof(small_sim_versions) / sizeof(small_sim_versions[0]); j++) {
            init_test(small_edge_file, small_demand_files[i], small_sim_versions[j]);

            for (int crossover_type = 0; crossover_type < CROSSOVER_TYPE; crossover_type++)
            {
                for (int mutation_type = 0; mutation_type < MUTATION_TYPE; mutation_type++)
                {
                    char filename[100];
                    sprintf(filename, "%s/outputs_sim%d_demand_%d_%d_%d.txt", folder_name, j + 1, i + 1, crossover_type, mutation_type);
                    freopen(filename, "a", stdout); // Redirect stdout to the file
                    printf("Crossover Type: %d, Mutation Type: %d\n", crossover_type, mutation_type);
                    // Call the test function to run the simulation
                    test(crossover_type, mutation_type);
                    fclose(stdout); // Close the file after writing
                }
            }
        }
    }

    return 0;
}