////////////////////////////////////////////////////////////////////////////////
//
// todoku : a general sudoku solver with an eternal TODO list...
//
// Copyright (c) Brett Hale 2019.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////


#include <stdexcept>
#include <string>
#include <cstdio>
#include <chrono>

#include <utility>
#include <bitset>
#include <memory>
#include <vector>
#include <list>
#include <algorithm>

#include <fstream>
#include <iostream> // (std::ifstream or std::cin)


using std::fprintf;


////////////////////////////////////////////////////////////////////////////////


class Game
{
public:

    static constexpr unsigned int max_size = (6);

    explicit Game (unsigned int n = 3);


    typedef std::pair<unsigned int, unsigned int> Address;

    static unsigned int row (const Address & address) {
        return std::get<0>(address); }
    static unsigned int col (const Address & address) {
        return std::get<1>(address); }

    class Cell
    {
    public:

        Cell (unsigned int n);

        Address address; // {{1, 1}, .., {n * n, n * n}}

        unsigned int potentials, solution;
        std::bitset<max_size * max_size> pset;


        void set (unsigned int); // {1, .., n * n}

        unsigned int row () const {
            return Game::row(address); } // {1, .., n * n}
        unsigned int col () const {
            return Game::col(address); } // {1, .., n * n}

        void dump () const;
    };


    unsigned int size; // (1 <= size (n) <= max_size)

    typedef std::vector<std::shared_ptr<Cell>> Group;
    std::vector<Group> grid; // row group state vector.


    bool solved () const;

    Cell & operator () (const Address & address) {
        return *(grid[row(address) - 1][col(address) - 1]); }

    const Cell & operator () (const Address & address) const {
        return *(grid[row(address) - 1][col(address) - 1]); }

    void dump () const;

};


////////////////////////////////////////////////////////////////////////////////


Game::Cell::Cell (unsigned int n)

    : potentials {n * n}, solution {n * n}
{
#if (0) // range check (1 <= n <= Game::max_size) :

    if (n < 1 || n > max_size)
        throw std::logic_error {"Cell ctor out of range\n"};
#endif

    if (max_size < (8))
        pset = (1ULL << potentials) - 1;
    else
    {
        for (unsigned int i = 0; i < potentials; i++)
            pset[i] = true;
    }
}


void
Game::Cell::set (unsigned int value)
{
    // TODO: what's the policy going to be here?

    potentials = 1, solution = value;

    if (max_size < (8))
        pset = 1ULL << (solution - 1);
    else
    {
        pset.reset(); pset[solution - 1] = true;
    }
}


void
Game::Cell::dump () const
{
    fprintf(stdout, "(%u, %u) ", row(), col());

    if (potentials == 1)
        fprintf(stdout, "solved: %u", solution);
    else
    {
        unsigned int i, nxn = solution;
        for (i = 0; i < nxn && !pset[i]; i++)
            ;
#if (0)
        if (i == nxn) // how did we get here?
            throw std::logic_error {"null potential set\n"};
#endif
        fprintf(stdout, ": %u", i + 1);
        for (i = i + 1; i < nxn; i++)
            if (pset[i]) fprintf(stdout, ", %u", i + 1);
    }

    fprintf(stdout, "\n");
}


Game::Game (unsigned int n)

    : size {n}, grid {(n * n), Group {(n * n), nullptr}}
{
#if (0) // range check (1 <= n <= Game::max_size) :

    if (n < 1 || n > max_size)
        throw std::logic_error {"Game ctor out of range\n"};
#endif

    unsigned int nxn = size * size;
    for (unsigned int i = 0; i < nxn; i++) // cell row:
    {
        for (unsigned int j = 0; j < nxn; j++) // cell col:
        {
            grid[i][j] = std::make_shared<Cell>(size);
            grid[i][j]->address = {i + 1, j + 1};
        }
    }
}


void
Game::dump () const
{
    unsigned int nxn = size * size;
    for (unsigned int i = 0; i < nxn; i++) // cell row:
    {
        for (unsigned int j = 0; j < nxn; j++) // cell col:
            grid[i][j]->dump();
    }
}


////////////////////////////////////////////////////////////////////////////////


bool
Game::solved () const
{
    decltype(Game::Cell::pset) group, none, all;
    unsigned int i, j, k, nxn = size * size;

    if (max_size < (8))
        all = (1ULL << nxn) - 1;
    else
        for (i = 0; i < nxn; i++) all[i] = true;

    for (i = 0; i < nxn; i++) // row (i) :
    {
        for (group = none, j = 0; j < nxn; j++)
        {
            const auto & cell = grid[i][j];
            if (cell->potentials != 1)
                return false;

            group[cell->solution - 1] = true;
        }

        if (group != all) return false;
    }

    for (j = 0; j < nxn; j++) // col (j) :
    {
        for (group = none, i = 0; i < nxn; i++)
        {
            const auto & cell = grid[i][j];
            group[cell->solution - 1] = true;
        }

        if (group != all) return false;
    }

    for (k = 0; k < nxn; k++) // box (k) :
    {
        unsigned int k_i = (k / size) * size;
        unsigned int k_j = (k % size) * size;

        for (group = none, i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                const auto & cell = grid[k_i + i][k_j + j];
                group[cell->solution - 1] = true;
            }
        }

        if (group != all) return false;
    }

    return true;
}


////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


static int
reduce (Game::Cell & cell, unsigned int value)
{
    if (cell.potentials == 1) // solved cell:
        return (cell.solution == value) ? (-1) : (0);

    if (!cell.pset[value - 1])
        return (0);

    cell.potentials--, cell.pset[value - 1] = false;

    if (cell.potentials == 1) // transition to solved cell:
    {
        unsigned int i, nxn = cell.solution;
        for (i = 0; i < nxn && !cell.pset[i]; i++)
            ;
#if (0)
        if (i == nxn) // how did we get here?
            throw std::logic_error {"null potential set\n"};
#endif
        cell.solution = i + 1;
    }

    return (+1); // i.e., potential set reduced.
}


static int
forward_solved (Game & game, Game::Cell & cell)
{
    int ret = 0;

    if (cell.potentials != 1) // unsolved cell:
        return (0);

    constexpr unsigned int max_peers =
        (3 * Game::max_size - 2) * Game::max_size - 1;

    Game::Address paddr[max_peers];
    unsigned int peers = 0; // until (3.n^2 - 2.n - 1)

    unsigned int n = game.size, value = cell.solution;
    unsigned int i, j, nxn = n * n;

    for (i = cell.row(), j = 1; j <= nxn; j++)
        if (j != cell.col())
            paddr[peers++] = {i, j};

    for (j = cell.col(), i = 1; i <= nxn; i++)
        if (i != cell.row())
            paddr[peers++] = {i, j};

    unsigned int box_i = ((cell.row() - 1) / n) * n;
    unsigned int box_j = ((cell.col() - 1) / n) * n;

    for (i = 1; i <= n; i++)
    {
        if (box_i + i == cell.row()) continue;
        for (j = 1; j <= n; j++)
        {
            if (box_j + j == cell.col()) continue;
            paddr[peers++] = {box_i + i, box_j + j};
        }
    }

    int rret, descent[max_peers];

    for (i = 0; i < peers; i++)
    {
        if ((rret = reduce(game(paddr[i]), value)) < 0)
            return (-1);

        descent[i] = rret, ret |= rret;
    }

    for (i = 0; ret >= 0 && i < peers; i++)
    {
        if (descent[i])
            ret |= forward_solved(game, game(paddr[i]));
    }

    return ret;
}


////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


static int
hidden_singles (Game & game, Game::Group & g)
{
    int ret = 0;

    unsigned int g_index[Game::max_size * Game::max_size];
    unsigned int nxn = game.size * game.size;

    unsigned int g_count[Game::max_size * Game::max_size];
    std::fill_n(g_count, nxn, 0);

    // on entry, assert all solved cells have been recursively
    // forwarded. the solution of a solved cell, therefore, can
    // not exist in more than one cell's p-set.

    // construct a histogram of p-set values in group:

    for (unsigned int i = 0; i < nxn; i++)
    {
        if (g[i]->potentials == 1)
        {
            unsigned int index = g[i]->solution - 1;
            g_count[index]++, g_index[index] = i;
        }
        else
        {
            unsigned int index = 0;
            for (; index < nxn; index++)
            {
                if (g[i]->pset[index])
                    g_count[index]++, g_index[index] = i;
            }
        }
    }

    for (unsigned int index = 0; index < nxn; index++)
    {
        if (g_count[index] == 0) // invalid state:
            return (-1);

        if (g_count[index] == 1)
        {
            auto & cell = *g[g_index[index]];
            if (cell.potentials == 1)
                continue;

            // assigning a hidden single value to a cell may
            // indirectly reduce p-sets of other group cells.

            cell.set(index + 1); ret = (+1);
            if (forward_solved(game, cell) < 0)
                return (-1);

            // this doesn't affect existing singular values
            // or invalid 'occupancy' (missing value) states.
        }
    }

    // note: tail-recursion can detect more than one singular
    // value in a single cell; this is itself an invalid state.

    return (ret > 0) ? (ret | hidden_singles(game, g)) : ret;
}


////////////////////////////////////////////////////////////////////////////////


constexpr unsigned int nxn_pairs (unsigned int n)
{
    return ((n * n) * (n * n - 1) / 2); // (n^2 pair combos)
}

#if (0) // inspired by Cantor's pairing function:

static unsigned int
forward_pairing (unsigned int x, unsigned int y)
{
    if ((0)) // (x <= y || y == 0)
        throw std::logic_error {"{x, y} s.t. x > y > 0\n"};

    return ((x - 1) * (x - 2) / 2 + (y - 1));
}

static std::pair<unsigned int, unsigned int>
inverse_pairing (unsigned int index)
{
    double det = 8.0 * index + 1.0;
    unsigned int w = static_cast<unsigned int>
        (std::floor((std::sqrt(det) - 1.0) / 2.0));

    unsigned int t = w * (w + 1) / 2;
    return {w + 2, index - t + 1}; // {x, y} s.t. x > y > 0.
}

#endif


static int
naked_pairs (Game & game, Game::Group & g)
{
    int ret = 0;

    for (unsigned int x = 2; x <= g.size(); x++)
    {
        for (unsigned int y = 1; y < x; y++) // x > y > 0:
        {
            unsigned int count = 0; // {x, y} cells.
            Game::Address paddr[2];

            for (const auto & cell : g)
            {
                if (cell->potentials != 2)
                    continue;

                if (cell->pset[x - 1] && cell->pset[y - 1])
                {
                    if (count == 2) // invalid state:
                        return (-1);

                    paddr[count++] = cell->address;
                }
            }

            if (count < 2) // no naked {x, y} pair:
                continue;

            for (auto & cell : g)
            {
                auto addr = cell->address;
                if (addr == paddr[0] || addr == paddr[1])
                    continue;

                int cret = 0;

                if ((cret |= reduce(*cell, x)) < 0)
                    return (-1);
                if ((cret |= reduce(*cell, y)) < 0)
                    return (-1);

                if (cret == 0) // x, y not in cell p-set:
                    continue;

                if ((cret |= forward_solved(game, *cell)) < 0)
                    return (-1);

                ret |= cret;
            }

            // since no statistics are retained, no assertions
            // are invalidated by continuing to comb for {x, y}
            // pairs - or restarting.
        }
    }

    return ret;
}


static int
hidden_pairs (Game & game, Game::Group & g)
{
    int ret = 0;

    unsigned int g_count[Game::max_size * Game::max_size];
    std::fill_n(g_count, g.size(), 0);

    constexpr auto max_pairs = nxn_pairs(Game::max_size);
    auto pairs = nxn_pairs(game.size);

    unsigned int p_count[max_pairs]; // {x, y} counts.
    std::fill_n(p_count, pairs, 0);

    Game::Address p_inverse[max_pairs]; // inverse pairings.

    // construct a histogram of p-set values and pairs:

    for (const auto & cell : g)
    {
        if (cell->potentials == 1) // no pairing possible:
        {
            g_count[cell->solution - 1]++;
            continue;
        }

        for (unsigned int i = 0; i < g.size(); i++)
            if (cell->pset[i]) g_count[i]++;

        for (unsigned int x = 2; x <= g.size(); x++)
        {
            if (!cell->pset[x - 1]) continue;
            for (unsigned int y = 1; y < x; y++) // x > y > 0:
            {
                if (!cell->pset[y - 1]) continue;

                auto index = (x - 1) * (x - 2) / 2 + (y - 1);
                p_count[index]++, p_inverse[index] = {x, y};
            }
        }
    }

    for (unsigned int p = 0; p < pairs; p++)
    {
        if (p_count[p] != 2) // no {x, y} -> (p) pairs:
            continue;

        auto x = p_inverse[p].first, y = p_inverse[p].second;
        if (g_count[x - 1] != 2 || g_count[y - 1] != 2)
            continue;

        for (auto & cell : g) // hidden pair:
        {
            if (!cell->pset[x - 1] || !cell->pset[y - 1])
                continue;

            if (cell->potentials == 2) // p-set = {x, y} :
                continue;

            for (unsigned int i = 1; i <= g.size(); i++)
            {
                if (i == x || i == y)
                    continue;

                if ((ret |= reduce(*cell, i)) < 0)
                    return (-1);
            }
        }

        // if (ret != 0) break;

        // exposing a hidden pair may 'invalidate' subsequent
        // hidden pairs. such false-positives are indicative of
        // overlapping pairs, and will fail the occupancy tests
        // in 'reduce' (above) or 'hidden_singles' (restart).
    }

    return ret;
}


////////////////////////////////////////////////////////////////////////////////


static int
group_deductions (Game & game, Game::Group & g)
{
    int ret = 0;

    // check that every value: {1, .., n * n} is present in each
    // group p-set union. if a value is unique to a cell, it has
    // been solved for that value, and is propagated:

#if (0)

    decltype(Game::Cell::pset) group, all;
    unsigned int nxn = game.size * game.size;

    if (Game::max_size < (8))
        all = (1ULL << nxn) - 1;
    else
    {
        for (unsigned int i = 0; i < nxn; i++)
            all[i] = true;
    }

    for (unsigned int i = 0; i < nxn; i++)
    {
        const auto & cell = g[i];
        if (cell->potentials != 1)
            break;

        group[cell->solution - 1] = true;
    }

    if (group == all) // self-consistent 'solved' group:
        return (0);
#endif

    // note: 'naked singles' are simply solved cells.

    if ((ret |= hidden_singles(game, g)) < 0)
        return (-1);

    // naked and hidden pair strategies:

    if ((ret |= naked_pairs(game, g)) < 0)
        return (-1);

    if ((ret |= hidden_pairs(game, g)) < 0)
        return (-1);

    return (ret > 0) ? (ret | group_deductions(game, g)) : ret;
}


static int
group_deductions (Game & game)
{
    int ret = 0;

    // each cell is a member of (3) groups: row, col, and box.

    unsigned int n = game.size, nxn = n * n;
    Game::Group g {(nxn), nullptr};

    for (unsigned int index = 0; index < nxn; index++)
    {
        for (unsigned int j = 0; j < nxn; j++)
            g[j] = game.grid[index][j];

        if ((ret |= group_deductions(game, g)) < 0) // (row)
            return (-1);

        for (unsigned int i = 0; i < nxn; i++)
            g[i] = game.grid[i][index];

        if ((ret |= group_deductions(game, g)) < 0) // (col)
            return (-1);

        unsigned int box_i = (index / n) * n;
        unsigned int box_j = (index % n) * n;

        for (unsigned int i = 0; i < n; i++)
        {
            for (unsigned int j = 0; j < n; j++)
                g[i * n + j] = game.grid[box_i + i][box_j + j];
        }

        if ((ret |= group_deductions(game, g)) < 0) // (box)
            return (-1);
    }

    return (ret > 0) ? (ret | group_deductions(game)) : ret;
}


////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


static int chute_deductions
(
    Game & game, Game::Group & chute, unsigned int x)
{
    unsigned int n = game.size;

    // classify: (0) ((x) not in box-line p-set), else (x)
    unsigned int classify[Game::max_size * Game::max_size];

    for (unsigned int i = 0; i < n; i++) // line (i) :
    {
        for (unsigned int j = 0; j < n; j++) // box (j) :
        {
            unsigned int box_line = i * n + j, found = 0;
            for (unsigned int k = 0; !found && k < n; k++)
            {
                if (chute[box_line * n + k]->pset[x - 1])
                    found = x;
            }

            classify[box_line] = found;
        }
    }

    // testing and updating the classification values yields
    // an extremely efficient test for box-line p-set unions,
    // with respect to (x), reflecting the current state:

    int ret = 0, cret;

    for (unsigned int i = 0; i < n; i++) // line (i) :
    {
        for (unsigned int j = 0; j < n; j++) // box (j) :
        {
            unsigned int box_line = i * n + j;
            if ((x = classify[box_line]) == 0)
                continue;

            bool pointing = true;
            for (unsigned int k = 0; pointing && k < n; k++)
            {
                if (k == i) continue;
                pointing = (classify[k * n + j] == 0);
            }

            for (unsigned int k = 0; pointing && k < n; k++)
            {
                if (k == j) continue;
                if (classify[box_line = i * n + k] == 0)
                    continue;

                for (unsigned int c = 0; c < n; c++)
                {
                    auto & cell = *chute[box_line * n + c];
                    if ((cret = reduce(cell, x)) > 0)
                        cret |= forward_solved(game, cell);

                    if ((ret |= cret) < 0) // invalid state:
                        return (-1);
                }

                classify[box_line] = 0; // update the box-line.
            }

            bool claiming = !pointing;
            for (unsigned int k = 0; claiming && k < n; k++)
            {
                if (k == j) continue;
                claiming = (classify[i * n + k] == 0);
            }

            for (unsigned int k = 0; claiming && k < n; k++)
            {
                if (k == i) continue;
                if (classify[box_line = k * n + j] == 0)
                    continue;

                for (unsigned int c = 0; c < n; c++)
                {
                    auto & cell = *chute[box_line * n + c];
                    if ((cret = reduce(cell, x)) > 0)
                        cret |= forward_solved(game, cell);

                    if ((ret |= cret) < 0) // invalid state:
                        return (-1);
                }

                classify[box_line] = 0; // update the box-line.
            }
        }
    }

    return ret;
}


static int
chute_deductions (Game & game)
{
    int ret = 0;

    unsigned int n = game.size, nxn = n * n;
    Game::Group chute {(nxn * n), nullptr};

    for (unsigned int floor = 0; floor < n; floor++)
    {
        for (unsigned int row = 0; row < n; row++)
        {
            unsigned int f_row = floor * n + row;
            for (unsigned int col = 0; col < nxn; col++)
                chute[row * nxn + col] = game.grid[f_row][col];
        }

        for (unsigned int x = 1; x <= nxn; x++)
            if ((ret |= chute_deductions(game, chute, x)) < 0)
                return (-1);
    }

    for (unsigned int tower = 0; tower < n; tower++)
    {
        for (unsigned int col = 0; col < n; col++)
        {
            unsigned int t_col = tower * n + col;
            for (unsigned int row = 0; row < nxn; row++)
                chute[col * nxn + row] = game.grid[row][t_col];
        }

        for (unsigned int x = 1; x <= nxn; x++)
            if ((ret |= chute_deductions(game, chute, x)) < 0)
                return (-1);
    }

    // return (ret > 0) ? (ret | chute_deductions(game)) : ret;

    // chute (locked candidate) deductions are too expensive
    // for exhaustive tail-recursion. any deductions (progress)
    // will be applied exhaustively to groups before repeating
    // chute deductions again.

    return ret;
}


////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////


struct Cell_Sort // tweaked sort functor:
{
    const Game & g;

    Cell_Sort (const Game & game) : g {game} {}

    bool operator () // bijection (1) -> (~0U) :
    (const Game::Address & x, const Game::Address & y) const {
        return ((g(x).potentials - 2) < (g(y).potentials - 2));
    }
};


static int
solve (Game & game, unsigned int level = 0)
{
    // fprintf(stdout, "*** solve: depth = %u\n\n", level);
    unsigned int n = game.size, nxn = n * n;

    if (level == 0) // recursively forward solved cells:
    {
        for (unsigned int i = 1; i <= nxn; i++)
        {
            for (unsigned int j = 1; j <= nxn; j++)
            {
                if (forward_solved(game, game({i, j})) < 0)
                    return (-1);
            }
        }
    }

    if (game.solved()) // immediate solution:
        return (0);

    for (int ret = (+1); ret > 0; ) // deductive loop:
    {
        // group_deductions is exhaustive:
        if ((ret = group_deductions(game)) < 0)
            return (-1);

        if (ret > 0 && game.solved())
            return (+1);

        // chute_deductions is *not* exhaustive:
        if ((ret = chute_deductions(game)) < 0)
            return (-1);

        if (ret > 0 && game.solved())
            return (+1);
    }

    // depth-first recursive search with trial candidates:

    // sort cell addresses by increasing p-set size, from (2)
    // to (n * n), with p-set sizes of (1) (representing solved
    // cells) being placed at the end (Cell_Sort) :

    std::vector<Game::Address> addr (nxn * nxn);

    for (unsigned int i = 0; i < nxn; i++) // cell row:
    {
        for (unsigned int j = 0; j < nxn; j++) // cell col:
            addr[i * nxn + j] = {i + 1, j + 1};
    }

    auto first = addr.begin(), last = addr.end();
    std::sort(first, last, Cell_Sort {game});

    Game trial {n};

    for (const auto & index : addr) // order: {2, .., nxn, 1}
    {
        auto & cell = game(index);
        for (unsigned int p = 1; p <= nxn; p++)
        {
            if (cell.potentials == 1) // no more candidates:
                break;

            if (!cell.pset[p - 1]) // (p) not in p-set:
                continue;

            for (unsigned int i = 1; i <= nxn; i++)
                for (unsigned int j = 1; j <= nxn; j++)
                    trial({i, j}) = game({i, j});

            auto & trial_cell = trial(cell.address);
            trial_cell.set(p);

            if (forward_solved(trial, trial_cell) < 0)
            {
                reduce(cell, p);
                continue;
            }

            if (solve(trial, level + 1) < 0)
            {
                reduce(cell, p);
                continue;
            }

            for (unsigned int i = 1; i <= nxn; i++)
                for (unsigned int j = 1; j <= nxn; j++)
                    game({i, j}) = trial({i, j});

            return (1); // (recursive) trial values solved game.
        }
    }

    return (-1);
}


////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

// game entry encodings : N = 1 (1x1 entry) to N = 6 (36x36 entries) :


typedef std::vector<std::pair<int, unsigned int>> map_n_container;

static const map_n_container map_n_1 {
    {'.',  0}, {'0',  0}, {'1',  1}};

static const map_n_container map_n_2 {
    {'.',  0}, {'0',  0}, {'1',  1}, {'2',  2}, {'3',  3}, {'4',  4}};

static const map_n_container map_n_3 {
    {'.',  0}, {'0',  0}, {'1',  1}, {'2',  2}, {'3',  3}, {'4',  4},
    {'5',  5}, {'6',  6}, {'7',  7}, {'8',  8}, {'9',  9}};

static const map_n_container map_n_4 {
    {'.',  0}, {'0',  1}, {'1',  2}, {'2',  3}, {'3',  4}, {'4',  5},
    {'5',  6}, {'6',  7}, {'7',  8}, {'8',  9}, {'9', 10}, {'A', 11},
    {'B', 12}, {'C', 13}, {'D', 14}, {'E', 15}, {'F', 16}};

static const map_n_container map_n_5 {
    {'.',  0}, {'0',  0}, {'A',  1}, {'B',  2}, {'C',  3}, {'D',  4},
    {'E',  5}, {'F',  6}, {'G',  7}, {'H',  8}, {'I',  9}, {'J', 10},
    {'K', 11}, {'L', 12}, {'M', 13}, {'N', 14}, {'O', 15}, {'P', 16},
    {'Q', 17}, {'R', 18}, {'S', 19}, {'T', 20}, {'U', 21}, {'V', 22},
    {'W', 23}, {'X', 24}, {'Y', 25}};

static const map_n_container map_n_6 {
    {'.',  0}, {'0',  1}, {'1',  2}, {'2',  3}, {'3',  4}, {'4',  5},
    {'5',  6}, {'6',  7}, {'7',  8}, {'8',  9}, {'9', 10}, {'A', 11},
    {'B', 12}, {'C', 13}, {'D', 14}, {'E', 15}, {'F', 16}, {'G', 17},
    {'H', 18}, {'I', 19}, {'J', 20}, {'K', 21}, {'L', 22}, {'M', 23},
    {'N', 24}, {'O', 25}, {'P', 26}, {'Q', 27}, {'R', 28}, {'S', 29},
    {'T', 30}, {'U', 31}, {'V', 32}, {'W', 33}, {'X', 34}, {'Y', 35},
    {'Z', 36}};

static const map_n_container *map_n[] =
{
    & map_n_1, & map_n_2, & map_n_3, & map_n_4, & map_n_5, & map_n_6
};


static std::vector<unsigned int>
map_game_data (unsigned int n, const std::string & s)
{
    std::vector<unsigned int> v;
    unsigned int nxn = n * n;

    if (s.length() < (nxn * nxn)) // return zero size:
        return v;

    // strictly speaking, this has O(n^6) lookups; in practice,
    // it's too fast to matter, and generates far less code than
    // O(n^4.log(n)) std::map implementations:

    v.reserve(nxn * nxn);
    for (unsigned int i = 0; i < nxn * nxn; i++)
    {
        bool found = false;
        for (auto p = map_n[n - 1]->cbegin();
             !found && p != map_n[n - 1]->cend(); ++p)
        {
            if (s[i] == p->first)
                found = true, v.push_back(p->second);
        }
    }

    return v;
}


static std::list<std::vector<unsigned int>>
read_game_data (unsigned int n, const char *file_name = nullptr)
{
    std::list<std::vector<unsigned int>> game_data;

    unsigned int nxn = n * n;
    std::string line;

    if (n > (sizeof(map_n) / sizeof(map_n[0])))
    {
        fprintf(stderr, "no game encoding for N = %u\n", n);
        return game_data; // (no game data)
    }

    if (file_name == nullptr) // single game from stdin:
    {
        std::getline(std::cin, line);

        if (line.length() == 0) // a 'blank' game:
            line.resize((nxn * nxn), '.');

        std::vector<unsigned int> v = map_game_data(n, line);

        if (v.size() == (nxn * nxn))
            game_data.push_back(std::move(v));
        else
            fprintf(stderr, "invalid game encoding\n");

        return game_data;
    }

    std::ifstream data_file {file_name};

    if (!data_file) // bad() or fail() bit:
    {
        fprintf(stderr, "could not open game data file: "
                "\"%s\"\n", file_name);
        return game_data; // (no game data)
    }

    while (std::getline(data_file, line))
    {
        if (line.length() == 0) // no 'blank' games:
            continue;
        if (line[0] == '#') // comment line:
            continue;

        std::vector<unsigned int> v = map_game_data(n, line);

        if (v.size() == (nxn * nxn))
            game_data.push_back(std::move(v));
    }

    if (game_data.empty())
        fprintf(stderr, "no valid games found for N = %u\n", n);

    return game_data;
}


////////////////////////////////////////////////////////////////////////////////


static int
print_grid (const Game & game)
{
    constexpr unsigned int max_side =
        (Game::max_size * (Game::max_size + 1) + 1);

    unsigned int n = game.size;
    unsigned int side = (n * (n + 1) + 1);

    char grid[max_side * (max_side + 1) + 1];
    unsigned int stride = side + 1;

    for (unsigned int i = 0; i < side; i++)
    {
        if (i % (n + 1) == 0)
        {
            for (unsigned int j = 0; j < side; j++)
                grid[i * stride + j] = ('-');
            for (unsigned int j = 0; j < side; j += n + 1)
                grid[i * stride + j] = ('+');
        }
        else
        {
            for (unsigned int j = 0; j < side; j += n + 1)
                grid[i * stride + j] = ('|');
        }

        grid[i * stride + side] = ('\n');
    }

    for (unsigned int i = 0; i < n * n; i++)
    {
        for (unsigned int j = 0; j < n * n; j++)
        {
            const auto & cell = game({i + 1, j + 1});
            unsigned int value = (cell.potentials == 1) ?
                cell.solution : (0);

            char c = '#'; // int -> char map fail.

            for (auto p = map_n[n - 1]->cbegin() + value;
                 c == '#' && p != map_n[n - 1]->cend(); ++p)
            {
                if (value == p->second)
                    c = static_cast<char>(p->first);
            }

            unsigned int row = (i + (i / n) + 1) * stride;
            grid[row + (j + (j / n) + 1)] = c;
        }
    }

    unsigned int count = side * stride;
    grid[count++] = ('\n');

    if (std::fwrite(grid, 1, count, stdout) != count)
        return (-1); // (write error)

    return (0);
}


////////////////////////////////////////////////////////////////////////////////


int main (int argc, char **argv)
{
    unsigned int n = (3); // the default sudoku size.

    if (argc > 1)
    {
        unsigned long u = 0;

        try
        {
            std::size_t nul_pos;
            u = std::stoul(argv[1], & nul_pos, 0);

            if (argv[1][nul_pos] != 0) // trailing...
                u = 0;
        }

        catch (const std::logic_error &)
        {
            u = 0;
        }

        if ((u < 1) || (u > Game::max_size))
        {
            fprintf(stderr, "todoku [size [file]], where: "
                    "size = 1 .. %u [%u]\n", Game::max_size, n);
            return (1);
        }

        n = static_cast<unsigned int>(u);
    }

    try
    {
        auto game_data =
            read_game_data(n, (argc > 2) ? argv[2] : nullptr);

        if (game_data.empty()) // reasons already given...
            return (1);

        typedef std::chrono::steady_clock timer;
        auto start = timer::now();

        for (const auto & v : game_data)
        {
            Game game {n};

            unsigned int nxn = n * n, clues = 0;
            for (unsigned int i = 0; i < nxn; i++)
            {
                for (unsigned int j = 0, u; j < nxn; j++)
                {
                    if ((u = v[i * nxn + j]) != 0)
                        clues++, game({i + 1, j + 1}).set(u);
                }
            }

            fprintf(stdout, "clues: %u\n\n", clues);
            print_grid(game);

            if (solve(game) < 0)
                fprintf(stdout, "invalid game\n\n");
            else
                print_grid(game);
        }

        auto finish = timer::now();
        std::chrono::duration<double> elapsed = finish - start;

        auto seconds = elapsed.count();
        fprintf(stdout, "elapsed time: %f seconds\n", seconds);

        auto games = static_cast<std::size_t>(game_data.size());
        fprintf(stdout, "games: %zu\n", games);
        auto average = seconds / static_cast<double>(games);
        fprintf(stdout, "average time: %f seconds\n", average);

        return (0);
    }

    catch (const std::logic_error & ex)
    {
        fprintf(stderr, "std::logic_error:\n");
        fprintf(stderr, "%s\n", ex.what());
    }

    catch (...)
    {
        fprintf(stderr, "exception caught (unspecified exception)\n");
    }

    return (1);
}


////////////////////////////////////////////////////////////////////////////////
