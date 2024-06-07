//
// Created by soumyadeep on 5/31/2024.
//

#ifndef NEWICKCOMPARATOR_NEWICK_COMPARATOR_H
#define NEWICKCOMPARATOR_NEWICK_COMPARATOR_H


#include <iostream>
#include <stdio.h>

using namespace std;

int CharToInt(char* c)
{
    int temp = 0;
    while (*c <= '9' && *c >= '0')
    {
        temp *= 10;
        temp += (*c - '0');
        scanf_s("%c", c, 1);
    }

    return temp;
}

template<class T>
struct ListElement//list element templa
{
    ListElement* next;
    T value;

    ListElement(T t)
    {
        value = t;
        next = NULL;
    }
};

template<class T>
struct List//lista template
{
    ListElement<T>* first;
    ListElement<T>* last;
    int size;
    List()
    {
        size = 0;
        first = NULL;
        last = NULL;
    }

    void Push(T c)
    {
        ListElement<T>* temp;
        temp = new ListElement<T>(c);

        if (size == 0)
        {
            first = temp;
            last = temp;
        }
        else
        {
            last->next = temp;
            last = temp;
        }

        size++;
    }

    T Get()
    {
        if (size > 0)return first->value;
        else return NULL;
    }

    T Get(int t)
    {
        ListElement<T>* temp = first;
        for (int i = 0; i < t; i++)
        {
            temp = temp->next;
        }

        return temp->value;
    }

    void Pop()
    {
        if (size > 0)
        {
            ListElement<T>* temp;
            temp = first;
            first = temp->next;
            size--;
            delete temp;
        }
    }

    ~List()
    {
        while (size > 0)
            Pop();
    }
};

struct Vertice;

struct Node
{
    int value;
    Node* brother;
    Node* son;
    Node* parent;
    Vertice* vertice;
    List<Node*>* list;

    Node(int val)
    {
        value = val;
        brother = NULL;
        son = NULL;
        parent = NULL;
        list = NULL;
        vertice = NULL;
    }

    Node(Node& v)
    {
        brother = v.brother;
        parent = v.parent;
        son = v.son;
        value = v.value;
        list = v.list;
        vertice = v.vertice;
    }

    void AddSon(Node* kid)
    {
        if (son == NULL) son = kid;
        else
        {
            if (brother == NULL)
            {
                Node* node = new Node(value);
                node->parent = parent;
                node->son = kid;
                brother = node;
            }
            else
            {
                if (brother->value == value)
                {
                    brother->AddSon(kid);
                }
                else
                {
                    Node* node = new Node(*this);
                    brother = node;
                }
            }
        }
    }

    void AddBrother(Node* bro)
    {
        if (brother == NULL)brother = bro;
        else brother->AddBrother(bro);
    }
};

struct Vertice
{
    int children;
    int name;
    List<Node*> nodes;
    List<Vertice*> vertices;
    Vertice* parent;
    int size;
    int* indeks;

    ~Vertice()
    {
        delete[] indeks;
    }

    Vertice()
    {
        children = 0;
        name = 0;
        size = 0;
        indeks = NULL;
        parent = NULL;
    }

    Vertice(char a)
    {
        name = a;
        children = 0;
        size = 0;
        indeks = NULL;
        parent = NULL;
    }

    void standardization()
    {
        size = nodes.size + vertices.size;
        indeks = new int[size];

        for (int i = 0; i < size; i++)
        {
            if (i < size - vertices.size)
            {
                indeks[i] = nodes.Get()->value;
                nodes.Pop();
            }
            else
            {
                indeks[i] = vertices.Get()->name;
                vertices.Pop();
            }
        }
    }
};


struct Tree
{
    Node* tip;
    int size;
    int ver;
    Vertice** inside;
    int** tab;
    List<Node*> nod;

    void Sort(List<Vertice*>* v)
    {
        int siz = ver;
        inside = new Vertice * [siz];
        v->Pop();
        for (int i = 0; i < siz; i++)
        {
            inside[i] = v->Get();
            v->Pop();
        }

        for (int i = 0; i < siz - 1; i++)
        {
            for (int j = 0; j < siz - 1 - i; j++)
            {
                if (inside[j]->children > inside[j + 1]->children) swap(inside[j], inside[j + 1]);
            }
        }
    }

    Tree()
    {
        tip = new Node(0);
        size = 0;
        ver = 0;
        inside = NULL;
        tab = NULL;
    }

    int fixVertice(Vertice* actual, List<Vertice*>* v)
    {
        v->Push(actual);
        for (int i = 0; i < actual->vertices.size; i++)
        {
            actual->children += fixVertice(actual->vertices.Get(i), v);
        }

        actual->children += actual->nodes.size;
        return actual->children;
    }

    void CreateTab()//tworzenie tablicy opisujacej ktory wierzcholek zawiera jakie wezly
    {
        tab = new int* [size];
        for (int i = 0; i < size; i++)
            tab[i] = new int[ver];

        for (int i = 0; i < ver; i++)
            for (int j = 0; j < size; j++)
                tab[j][i] = 0;

        for (int i = 0; i < ver; i++)
        {
            for (int j = 0; j < inside[i]->size; j++)
            {
                if (inside[i]->indeks[j] > 0) tab[inside[i]->indeks[j] - 1][i] = 1;
                else
                {
                    for (int k = 0; k < size; k++)
                    {
                        if (tab[k][(inside[i]->indeks[j] + 1) * -1] == 1) tab[k][i] = 1;
                    }
                }
            }
        }
    }

    void CreateTree()
    {
        char c = 0;
        Node* actual = tip;
        Vertice* actualver;
        actualver = new Vertice();
        scanf_s("%c", &c, 1);
        while (c != ';')
        {
            if (c <= '9' && c >= '0')
            {
                size++;
                if (actual->value != 0)
                {
                    Node* temp;
                    temp = new Node(CharToInt(&c));
                    nod.Push(temp);
                    actual->AddBrother(temp);
                    temp->parent = actual->parent;
                    temp->vertice = actual->vertice;
                    actual = temp;
                }
                else actual->value = CharToInt(&c);

                actualver->nodes.Push(actual);
            }
            else
            {
                if (c == '(')
                {
                    Node* temp = new Node(0);
                    nod.Push(temp);
                    actual->AddSon(temp);
                    temp->parent = actual;
                    temp->vertice = new Vertice('a' + ver);
                    actualver->vertices.Push(temp->vertice);
                    temp->vertice->parent = actualver;
                    actualver = temp->vertice;
                    actual = temp;
                    ver++;
                }
                else if (c == ')')
                {
                    actual = actual->parent;
                    actualver = actualver->parent;
                }

                scanf_s("%c", &c, 1);
            }
        }
        List<Vertice*> v;
        actualver->children = fixVertice(actualver, &v);

        Sort(&v);

        for (int i = 0; i < ver; i++)
        {
            inside[i]->standardization();
            inside[i]->name = -(i + 1);
        }

        CreateTab();
    }

    ~Tree()
    {
        while (nod.size > 0)
        {
            Node* n;
            n = nod.Get();
            nod.Pop();
            delete n;
        }

        for (int i = 0; i < ver; i++)
            delete inside[i];

        delete[] tab;
        delete inside;
        delete tip;
    }

};

void func(int** search, int width, int height, bool rotate, int* tab, int i, int j, int sum, int* maks)
{
    if (j < height)
    {
        for (int k = 0; k < width; k++)
        {
            bool cos = true;
            for (int g = 0; g < j; g++)
                if (tab[g] == k) cos = false;

            if (cos)
            {
                int temp = 0;
                if (rotate) temp = search[k][j];
                else temp = search[j][k];
                int* tabb = new int[height];
                for (int h = 0; h < height; h++)
                    tabb[h] = tab[h];

                tabb[j] = k;
                func(search, width, height, rotate, tabb, k, j + 1, sum + temp, maks);
            }
        }
    }
    else
    {
        if (sum > *maks) *maks = sum;

        delete[] tab;
    }
}


int MaksLine(int** search, int width, int height, bool rotate)
{
    int maks = 0;
    int* tab;
    tab = new int[height];

    for (int i = 0; i < width; i++)
    {
        tab[0] = i;
        int sum;
        if (rotate) sum = search[i][0];
        else sum = search[0][i];

        func(search, width, height, rotate, tab, 0, 1, sum, &maks);
    }

    delete[] tab;
    return maks;
}


int isomorfic(Vertice* first, Vertice* second, int** tab, int** ftab, int** stab)
{
    //najlepiej te taclice zrobic raz dla kazdego drzewa wiec do konstruktora trzeba to wrzucic pozniejk
    int maks = 0;
    int** search;
    search = new int* [first->size];
    for (int i = 0; i < first->size; i++)
        search[i] = new int[second->size];

    for (int i = 0; i < first->size; i++)
        for (int j = 0; j < second->size; j++)
            search[i][j] = 0;

    for (int i = 0; i < first->size; i++)
    {
        for (int j = 0; j < second->size; j++)
        {
            if (first->indeks[i] > 0 && second->indeks[j] > 0)
            {
                if (first->indeks[i] == second->indeks[j])search[i][j] = 1;
                else search[i][j] = 0;
            }
            else
            {
                if (first->indeks[i] < 0 && second->indeks[j] < 0) search[i][j] = tab[(first->indeks[i] + 1) * -1][(second->indeks[j] + 1) * -1];
                else
                {
                    if (first->indeks[i] > 0)
                    {
                        if (stab[first->indeks[i] - 1][(second->indeks[j] + 1) * -1] == 1) search[i][j] = 1;
                        else search[i][j] = 0;
                    }
                    else
                    {
                        if (ftab[second->indeks[j] - 1][(first->indeks[i] + 1) * -1] == 1) search[i][j] = 1;
                        else search[i][j] = 0;
                    }
                }
            }
        }
    }
    //wyliczanie maksymalnej liczby z tabeli
    if (first->size > second->size) maks = MaksLine(search, first->size, second->size, true);
    else maks = MaksLine(search, second->size, first->size, false);

    //sprawdzanie dzieci i porownywanie ich z aktualnym wezlem
    for (int i = 0; i < first->size; i++)
    {
        if (first->indeks[i] < 0)
        {
            if (maks < tab[(first->indeks[i] + 1) * -1][(second->name + 1) * -1])maks = tab[(first->indeks[i] + 1) * -1][(second->name + 1) * -1];
        }
    }

    for (int i = 0; i < second->size; i++)
    {
        if (second->indeks[i] < 0)
        {
            if (maks < tab[(first->name + 1) * -1][(second->indeks[i] + 1) * -1])maks = tab[(first->name + 1) * -1][(second->indeks[i] + 1) * -1];
        }
    }

    delete[] search;
    return maks;
}

int Calculate(Tree* one, Tree* two)
{
    int size = one->size;
    int maks = 0;
    int** tab;
    tab = new int* [one->ver];
    for (int i = 0; i < one->ver; i++)
        tab[i] = new int[two->ver];

    for (int i = 0; i < one->ver; i++)
        for (int j = 0; j < two->ver; j++)
            tab[i][j] = 0;

    for (int i = 0; i < one->ver + two->ver; i++)
    {
        for (int j = 0; j <= i && j < one->ver; j++)
        {
            for (int k = 0; k <= i && k < two->ver; k++)
            {
                if (j + k == i)
                {
                    tab[j][k] = isomorfic(one->inside[j], two->inside[k], tab, one->tab, two->tab);
                    if (tab[j][k] > maks)maks = tab[j][k];
                }
            }
        }
    }

    delete[] tab;
    //cout << size - maks << endl;
    return size - maks;
}

#endif //NEWICKCOMPARATOR_NEWICK_COMPARATOR_H
