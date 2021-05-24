#include <cstdlib>
#include <ctime>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <fstream>

long long GCD(long long a, long long b)
{
    while (b != 0)
    {
        long long r = a % b;
        a = b;
        b = r;
    }
    return a;
}
long long Mod(long long a, long long b)
{
    if (b < 0) return -Mod(-a, -b);
    long long res = a % b;
    return res < 0 ? res + b : res;
}

long long GCD_Ext(long long a, long long b, long long& x, long long& y)
{
    if (a == 0) { x = 0; y = 1; return b; }
    long long x1 = 0, y1 = 0;
    long long d = GCD_Ext(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}

long long find_x(long long* equations, long long N, long long& modulus)
{
    long long modulo = 1;
    for (long long i = 0; i < N; i++)
    {
        if (equations[i * 3] != 1)
        {
            if (equations[i * 3 + 1] >= equations[i * 3 + 2])
                equations[i * 3 + 1] = Mod(equations[i * 3 + 1], equations[i * 3 + 2]);
            long long inv1, inv2;
            GCD_Ext(equations[i * 3], equations[i * 3 + 2], inv1, inv2);
            equations[i * 3 + 1] = Mod(equations[i * 3 + 1] * Mod(inv1, equations[i * 3 + 2]), equations[i * 3 + 2]);
        }
        modulo *= equations[i * 3 + 2];
    }
    modulus = modulo;

    long long* inv = new long long[N];
    for (int i = 0; i < N; i++)
    {
        long long x, y;
        GCD_Ext(modulo / equations[i * 3 + 2], equations[i * 3 + 2], x, y);
        x = Mod(x, equations[i * 3 + 2]);
        inv[i] = x;
    }

    long long res = 0;

    for (int i = 0; i < N; i++)
        res += Mod(equations[i * 3 + 1] * inv[i] * (modulo / equations[i * 3 + 2]), modulo);
    return Mod(res, modulo);
}

bool check(long long* equations, long long N)
{
    for (long long i = 0; i < N-1; i++)
        for (long long j = i+1; j < N; j++)
            if (GCD(equations[3 * j + 2], equations[3 * i + 2]) != 1) return false;

    return true;
}


void main(int argc, char* argv[])
{

    int ProcNum, ProcRank;

    long long N, N_proc, modulus, x, x_, mod_;
    double time,time2;
    long long* equations; long long* eqs; long long* eq_f;
    const int one = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Status status;

    if (ProcRank == 0)
    {
        std::cout << "Enter n:\n";
        do
            std::cin >> N;
        while (N % ProcNum != 0 || N < ProcNum);
    }

    MPI_Bcast(&N, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    equations = (long long*)malloc((N * 3) * sizeof(long long));
    char y = NULL;
    if (ProcRank == 0)
    {
        std::cout << "File input or manual? f/m\n";
        std::cin >> y;
        if (y == 'm')
        {
            for (int i = 0; i < N; i++)
            {
                std::cout << "Enter " << i + 1 << "th equation coefficients:\n";
                std::cin >> equations[i * 3] >> equations[i * 3 + 1] >> equations[i * 3 + 2];
            }
        }
        else
            if (y == 'f')
            {
                std::ifstream in;
                std::string file;
                while (!in.is_open())
                {
                    std::cout << "Enter file name.\n";
                    std::cin >> file;
                    in.open("C:\\temp\\" + file);
                }

                for (int i = 0; i < N; i++)
                    {
                        long long temp1, temp2, temp3;
                        in >> temp1 >> temp2 >> temp3;
                        equations[3 * i] = temp1;
                        equations[3 * i + 1] = temp2;
                        equations[3 * i + 2] = temp3;
                        //in >> equations[3 * i] >> equations[3 * i + 1] >> equations[3 * i + 2];
                    }

                in.close();
            }
            else
            {
                std::cout << "Wrong input, program shall exit.\n";
                MPI_Abort(MPI_COMM_WORLD, -2);
            }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    N_proc = N / ProcNum;
    eqs = (long long*)malloc((N_proc * 3) * sizeof(long long));
    MPI_Scatter(equations, N_proc * 3, MPI_LONG_LONG, eqs, N_proc * 3, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    eq_f = (long long*)malloc((ProcNum * 3) * sizeof(long long));

    if (ProcRank == 0) time = MPI_Wtime();
    if (!check(eqs, N_proc))
        MPI_Abort(MPI_COMM_WORLD, -1);
    x = find_x(eqs, N_proc, modulus);
    if (ProcRank == 0) time = MPI_Wtime()-time;

    if (ProcRank != 0)
    {
        MPI_Send(&x, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&modulus, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    }if (ProcRank == 0)
    {
        eq_f[0] = 1;
        eq_f[1] = x;
        eq_f[2] = modulus;
        for (int i = 1; i < ProcNum; i++) {
            MPI_Recv(&x, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&modulus, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
            eq_f[3 * i] = 1;
            eq_f[3 * i + 1] = x;
            eq_f[3 * i + 2] = modulus;
        }
        time2 = MPI_Wtime();
        if (!check(equations, N))
             MPI_Abort(MPI_COMM_WORLD, -1);
        x_ = find_x(eq_f, ProcNum, mod_);
        time2 = MPI_Wtime() - time2;
        std::cout << "Res is " << x_ << " mod " << mod_ << " - calculated parallel for " << time+time2 << "\n";

    }
    delete[] equations;
    delete[] eqs;
    delete[] eq_f;
    MPI_Finalize();
}