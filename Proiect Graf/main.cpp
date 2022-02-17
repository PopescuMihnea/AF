#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <queue>
#include <climits>
#define PROBLEMA 2
#define inf INT_MAX-1000000
using namespace std;
ifstream f;
ofstream g;
class graf
{
    unsigned int n,m;
    vector< vector<pair<int,int> > > lista_vecini;
    void BFS(unsigned int,int*);
    void DFS(unsigned int, bool*,unsigned int&);
    void biconex(unsigned int, int*, int*, unsigned int&, stack<pair<unsigned int, unsigned int>>&, vector<unsigned int>&);
    void ctc(unsigned int, int*, int*, bool*, unsigned int&, unsigned int&, stack<unsigned int>&, vector <unsigned int>&);
    void sortaret(unsigned int, bool*, stack<unsigned int>&);
    vector<vector<int>> criticalConnections(int,vector<vector<int>>&);
    bool bfs_cuplaj(int*, int*, int*);
    bool dfs_cuplaj(unsigned int, int*, int*, int*);
public:
    graf(const unsigned int &);
    ~graf();
    void addMuchie(const unsigned int &, const unsigned int &, const int &);
    void setn(const unsigned int &);
    void setm(const unsigned int &);
    int getn();
    int getm();
    void getvecini();
    void sortvecini();
//de aici incepe numerotarea problemelor (de la 1)
    void StartBFS();
    void StartDFS();
    void Startbiconex();
    void Startctc();
    void Startsortaret();
    void havelhakimi();
    void StartcriticalConnections();
    void apm();
    void disjoint();
    void dijkstra();
    void Bellman_Ford();
    void maxflow();
    void royfloyd();
    void darb();
    void ciclueuler();
    void cicluhamilton();
    void cuplaj();
};

/*
0(gol):graf construit manual
fara cost
1:graf orientat cu nod de start
2:graf neorientat
3:graf orientat
6:graf neorientat unde se da doar nr de noduri
cu cost
4:graf orientat
5:graf neorientat
*/
//se pot seta manual nr de muchii,noduri cu seteri
//se pot adauga muchii cu metoda addMuchie
graf::graf(const unsigned int &opt=0)
{

    switch (opt)
    {
    case 1:
    case 3:
    case 4:
    {
        unsigned int x,y,cost;
        f>>this->n>>this->m;
        if (opt==1)
            f>>x;
        this->lista_vecini.resize(n+1);
        for (unsigned int i=0; i<this->m; ++i)
        {
            f>>x>>y;
            if (opt==1 || opt==3)
                this->lista_vecini[x].push_back(make_pair(y,0));
            else
            {
                f>>cost;
                this->lista_vecini[x].push_back(make_pair(y,cost));
            }
        }
        break;
    }
    case 2:
    case 5:
    case 6:
    {
        unsigned int x,y,cost;
        if (opt!=6)
            f>>this->n>>this->m;
        else
        {
            f>>this->n;
            this->m=this->n-1;
        }
        this->lista_vecini.resize(n+1);
        for (unsigned int i=0; i<this->m; ++i)
        {
            f>>x>>y;
            if (opt!=5)
            {
                this->lista_vecini[x].push_back(make_pair(y,0));
                this->lista_vecini[y].push_back(make_pair(x,0));
            }
            else
            {
                f>>cost;
                this->lista_vecini[x].push_back(make_pair(y,cost));
                this->lista_vecini[y].push_back(make_pair(x,cost));
            }

        }
        break;
    }
    default:
        break;
    }
}

graf::~graf()
{
    //simbolic
    this->n=this->m=0;
    for (size_t i=0; i<this->lista_vecini.size(); ++i)
        this->lista_vecini[i].clear();
    this->lista_vecini.clear();
}

void graf::setn(const unsigned int &n)
{
    this->n=n;
}

void graf::setm(const unsigned int &m)
{
    this->m=m;
}

//adauga muchie de la x la y cu costul c, costul este optional
void graf::addMuchie(const unsigned int &x, const unsigned int &y, const int &c=0)
{
    if (this->lista_vecini.size()!=this->n+1)
        this->lista_vecini.resize(this->n+1);
    this->lista_vecini[x].push_back(make_pair(y,c));
}

int graf::getm()
{
    return this->m;
}

int graf::getn()
{
    return this->n;
}

//afiseaza lista de vecini
void graf::getvecini()
{
    for (unsigned int i=1-min(1,int(lista_vecini[0].size())); i<this->lista_vecini.size()-min(1,int(lista_vecini[0].size())); ++i)
    {
        cout<<i<<": ";
        for (unsigned int j=0; j<this->lista_vecini[i].size(); ++j)
            cout<<this->lista_vecini[i][j].first<<" ";
        cout<<"\n";
    }
}

//afiseaza vecinii sortati lexicografic
void graf::sortvecini()
{
    unsigned int j,k;
    unsigned int ap[this->n+1];
    for (j=0; j<=this->n; ++j)
        ap[j]=0;
    for (unsigned int i=1; i<this->lista_vecini.size(); ++i)
    {
        k=0;
        for (j=0; j<this->lista_vecini[i].size(); ++j)
            ++ap[this->lista_vecini[i][j].first];
        for (j=1; j<=this->n; ++j)
        {
            if (ap[j])
            {
                this->lista_vecini[i][k++].first=j;
                ap[j]=0;
            }
        }
    }
}

void graf::StartBFS()
{
    ifstream h;
    unsigned int start;
    h.open("bfs.in",std::ifstream::in);
    h>>start>>start>>start;
    h.close();
    int d[this->n+1];
    for (unsigned int i=1; i<=this->n; ++i)
        d[i]=-1;
    BFS(start,d);
    for (unsigned int i=1; i<=this->n; ++i)
    {
        g<<d[i]<<" ";
    }
}

//afiseaza distanta de la un nod de start la celelalte noduri
//parcurgem vecinii nevizitati ai nodurilor dintr-o coada, incrementam distanta fiului cu 1 fata
//de cea a tatalui si adaugam fiul in coada
void graf::BFS(unsigned int nod_start,int d[])
{
    unsigned int index=0;
    vector <int> coada;
    coada.push_back(nod_start);
    d[nod_start]=0;
    while (index<coada.size())
    {
        int nod_curent=coada[index++];
        for (unsigned int i=0; i<this->lista_vecini[nod_curent].size(); ++i)
            if (d[this->lista_vecini[nod_curent][i].first]==-1)
            {
                coada.push_back(this->lista_vecini[nod_curent][i].first);
                d[this->lista_vecini[nod_curent][i].first]=d[nod_curent]+1;

            }
    }
}

void graf::StartDFS()
{
    bool viz[this->n+1];
    unsigned int nrconex;
    for (unsigned int i=1; i<=this->n; ++i)
        viz[i]=false;
    nrconex=0;
    for (unsigned int i=1; i<=this->n; ++i)
        if (!viz[i])
        {
            nrconex++;
            DFS(i,viz,nrconex);
        }
    g<<nrconex;

}

//gaseste nr de comp conexe
//aceiasi paradigma ca la bfs doar ca imediat ce gasim un vecin nevizitat
//il vizitam si il procesam, nodul precedent fiind prelucrat in continuare dupa ce se goleste stiva pana la el
void graf::DFS(unsigned int nod_curent,bool viz[],unsigned int &nrconex)
{
    viz[nod_curent]=true;
    for (unsigned int i=0; i<this->lista_vecini[nod_curent].size(); ++i)
        if (!viz[this->lista_vecini[nod_curent][i].first])
            DFS(this->lista_vecini[nod_curent][i].first,viz,nrconex);
}

void graf::Startbiconex()
{
    int niv[this->n+1];
    int niv_int[this->n+1];
    unsigned int nrbiconex=0;
    stack <pair<unsigned int,unsigned int>> stiva;
    vector <unsigned int> sol;
    for (unsigned int i=1; i<=this->n; ++i)
    {
        niv[i]=-1;
        niv_int[i]=0;
    }
    niv[1]=0;
    biconex(1,niv,niv_int,nrbiconex,stiva,sol);
    g<<nrbiconex<<"\n";
    for (auto i = sol.begin(); i != sol.end(); ++ i)
        if (*i!=0)
            g<<*i<<" ";
        else
            g<<"\n";
}

//arata numarul si muchiile componentelor biconexe
//construim un arbore al parcurgerii dfs cu drumuri de intoarcere (in sus) catre noduri deja
//vizitate.daca un nod are un vecin care are nivelul de int mai mare decat nivelul sau
//inseamna ca este nod critic (ca un cap de pod) si face parte din muchia critica, vecinul fiind celalalt capat
//cand gasim un nod critic dam pop la muchiile care sunt salvate intr-o stiva in ordinea parcurgerii pana cand
//dam de nodurile ce formeaza muchia critica (inclusiv ele)
void graf::biconex(unsigned int nod_curent,int niv[],int niv_int[],unsigned int &nrbiconex,stack <pair<unsigned int,unsigned int>> &stiva,vector <unsigned int> &sol)
{
    for (unsigned int i=0; i<this->lista_vecini[nod_curent].size(); ++i)
        if (niv[this->lista_vecini[nod_curent][i].first]==-1)
        {
            niv[this->lista_vecini[nod_curent][i].first]=niv[nod_curent]+1;
            niv_int[this->lista_vecini[nod_curent][i].first]=niv[nod_curent]+1;
            stiva.push(make_pair(nod_curent,this->lista_vecini[nod_curent][i].first));
            biconex(this->lista_vecini[nod_curent][i].first,niv,niv_int,nrbiconex,stiva,sol);
            niv_int[nod_curent]=min(niv_int[nod_curent],niv_int[this->lista_vecini[nod_curent][i].first]);
            if (niv_int[this->lista_vecini[nod_curent][i].first]>=niv[nod_curent])
            {
                nrbiconex++;
                unordered_map<unsigned int,bool> ap;
                unsigned int aux1;
                unsigned int aux2;
                do
                {

                    aux1=stiva.top().first;
                    aux2=stiva.top().second;
                    if (!ap[aux1])
                    {
                        sol.push_back(aux1);
                        ap[aux1]=1;
                    }
                    if (!ap[aux2])
                    {
                        sol.push_back(aux2);
                        ap[aux2]=1;
                    }
                    stiva.pop();
                }
                while(aux1!=nod_curent || aux2!=this->lista_vecini[nod_curent][i].first);
                sol.push_back(0);
            }
        }
        else if (niv[nod_curent]-1!=niv[this->lista_vecini[nod_curent][i].first])
            niv_int[nod_curent]=min(niv_int[nod_curent],niv[this->lista_vecini[nod_curent][i].first]);
    //scadem 1 ca sa nu luam in considerare si tatal nodului
}

void graf::Startctc()
{
    int indx[this->n+1];
    int indx_min[this->n+1];
    bool peStiva[this->n+1];
    unsigned int nrctc=0;
    unsigned int it=1;
    unsigned int i;
    stack <unsigned int> stiva;
    vector <unsigned int> sol;
    for (i=1; i<=this->n; ++i)
    {
        indx[i]=-1;
        indx_min[i]=0;
        peStiva[i]=false;
    }
    for(i=1; i<=this->n; ++i)
        if (indx[i]==-1)
            ctc(i,indx,indx_min,peStiva,nrctc,it,stiva,sol);
    g<<nrctc<<"\n";
    for (auto i = sol.begin(); i != sol.end(); ++ i)
        if (*i!=0)
            g<<*i<<" ";
        else
            g<<"\n";

}

//arata numarul si nodurile componentelor tare conexe
//aceiasi tehnica de la biconex, doar ca luam indexul nodului in loc de nivel
//daca nodul are index de intoarcere egal cu indexul sau inseamna ca este radacina comp conexe
//asadar golim stiva de toate nodurile pana dam de radacina (inclusiv radacina)
//daca un nod nu e pe stiva inseamna ca deja e intr-o ctc si il ignoram (peStiva e ca un fel de viz)
void graf::ctc(unsigned int nod_curent,int indx[], int indx_min[], bool peStiva[], unsigned int &nrctc, unsigned int &it, stack<unsigned int> &stiva, vector<unsigned int> &sol)
{
    indx[nod_curent]=it;
    indx_min[nod_curent]=it++;
    stiva.push(nod_curent);
    peStiva[nod_curent]=true;
    for (unsigned int i=0; i<this->lista_vecini[nod_curent].size(); ++i)
        if (indx[this->lista_vecini[nod_curent][i].first]==-1)
        {
            ctc(this->lista_vecini[nod_curent][i].first,indx,indx_min,peStiva,nrctc,it,stiva,sol);
            indx_min[nod_curent]=min(indx_min[nod_curent],indx_min[this->lista_vecini[nod_curent][i].first]);
        }
        else if (peStiva[this->lista_vecini[nod_curent][i].first])
            indx_min[nod_curent]=min(indx_min[nod_curent],indx[this->lista_vecini[nod_curent][i].first]);
    if (indx_min[nod_curent]==indx[nod_curent])
    {
        ++nrctc;
        unsigned int nod_aux;
        do
        {
            nod_aux=stiva.top();
            peStiva[nod_aux]=false;
            sol.push_back(nod_aux);
            stiva.pop();
        }
        while (nod_aux!=nod_curent);
        sol.push_back(0);
    }
}

void graf::Startsortaret()
{
    bool viz[this->n+1];
    stack <unsigned int> sol;
    unsigned int i;
    for (i=1; i<=this->n; ++i)
        viz[i]=false;
    for (i=1; i<=this->n; ++i)
        if (!viz[i])
            sortaret(i, viz, sol);
    while (!sol.empty())
    {
        //afisam in postordine invers
        g<<sol.top()<<" ";
        sol.pop();
    }
}

//sorteaza topologic nodurile unui graf
//trucul e sa afisam nodurile in postordinea unei parcurgeri DFS in ordine inversa
void graf::sortaret(unsigned int nod_curent, bool viz[], stack<unsigned int> &sol)
{
    viz[nod_curent]=true;
    for (unsigned int i=0; i<this->lista_vecini[nod_curent].size(); ++i)
        if (!viz[lista_vecini[nod_curent][i].first])
        {
            sortaret(lista_vecini[nod_curent][i].first, viz, sol);
        }
    //introducem nodurile in postordine (dupa ce ies din stiva) in stiva solutie
    sol.push(nod_curent);
}

//arata daca exista un graf neorientat cu gradele date
//sortam descrescator dupa grade
//eliminam gradul cel mai mare(grad_max) si din urmatoarele grad_max grade scadem 1
//daca dam de o valoare negativa sau nu mai avem din ce sa scadem nu se poate forma graful
void graf::havelhakimi()
{
    bool zero=false;
    int vmax=0;
    vector <int> grad;
    unsigned int x;
    while (f>>x)
        grad.push_back(x);
    while (!zero)
    {
        zero=true;
        int ap[grad.size()+1];
        int k=0;
        for (unsigned int i=0; i<=grad.size(); ++i)
            ap[i]=0;
        for (unsigned int i=0; i<grad.size(); ++i)
        {
            ++ap[grad[i]];
            vmax=max(vmax,grad[i]);
        }
        for (unsigned int i=0; i<=vmax; ++i)
            for (unsigned int j=0; j<ap[i]; ++j)
                grad[k++]=i;
        unsigned int nrs=grad[grad.size()-1];
        grad.pop_back();
        if (nrs>grad.size())
        {
            cout<<"Nu";
            g<<"Nu";
            return;
        }
        for (unsigned int i=grad.size()-1; (int)(i)>(int)(grad.size()-nrs-1); --i)
        {
            --grad[i];
            if (grad[i]<0)
            {
                cout<<"Nu";
                g<<"Nu";
                return;
            }
            if (grad[i]!=0)
                zero=false;
        }
    }
    cout<<"Da";
    g<<"Da";
}

void graf::StartcriticalConnections()
{
    vector<vector<int>> muchii;
    int n,x,y,m;
    unsigned int nrt;
    f>>nrt;
    for (unsigned int t=0; t<nrt; ++t)
    {
        muchii.resize(0);
        f>>n>>m;
        for (int i=0; i<m; ++i)
        {
            f>>x>>y;
            muchii.push_back(vector<int> {x,y});
        }
        vector<vector<int>> sol(this->criticalConnections(n,muchii));
        for (size_t i=0; i<sol.size(); ++i)
        {
            cout<<"("<<sol[i][0]<<","<<sol[i][1]<<") ";
            g<<"("<<sol[i][0]<<","<<sol[i][1]<<") ";
        }
        cout<<"\n\n";
        g<<"\n\n";
    }
}

//arata muchiile critice
//ca la biconex doar ca aratam doar muchiile critice
vector<vector<int>> graf::criticalConnections(int n, vector<vector<int>>& connections)
{
    unsigned int i;
    this->setn(n);
    this->setm(connections.size());
    stack <unsigned int> stiva;
    vector<vector<int>> sol;
    unsigned int icopil[n];
    this->lista_vecini.clear();
    this->lista_vecini.resize(n);
    int niv[n];
    int niv_int[n];
    for (int i=0; i<n; ++i)
    {
        niv[i]=-1;
        niv_int[i]=0;
        icopil[i]=0;
    }
    for (i=0; i<connections.size(); ++i)
    {
        this->lista_vecini[connections[i][0]].push_back(make_pair(connections[i][1],0));
        this->lista_vecini[connections[i][1]].push_back(make_pair(connections[i][0],0));
    }
    niv[0]=1;
    stiva.push(0);
    while (!stiva.empty())
    {
        unsigned int nod_curent=stiva.top();
        for (i=icopil[nod_curent]; i<this->lista_vecini[nod_curent].size(); ++i)
        {
            ++icopil[nod_curent];
            if (niv[this->lista_vecini[nod_curent][i].first]==-1)
            {
                niv[this->lista_vecini[nod_curent][i].first]=niv[nod_curent]+1;
                niv_int[this->lista_vecini[nod_curent][i].first]=niv[nod_curent]+1;
                stiva.push(this->lista_vecini[nod_curent][i].first);
                nod_curent=this->lista_vecini[nod_curent][i].first;
                break;
            }
            else if (niv[nod_curent]-1!=niv[this->lista_vecini[nod_curent][i].first])
                niv_int[nod_curent]=min(niv_int[nod_curent],niv[this->lista_vecini[nod_curent][i].first]);
        }
        if (icopil[nod_curent]==this->lista_vecini[nod_curent].size())
        {
            unsigned int copil=stiva.top();
            stiva.pop();
            if (!stiva.empty())
            {
                niv_int[stiva.top()]=min(niv_int[copil],niv_int[stiva.top()]);
                if (niv_int[copil]>niv[stiva.top()])
                    sol.push_back({int(stiva.top()),int(copil)});
            }

        }
    }
    return sol;
}

//construieste si arata costul unui arbore partial de cost minim
//e folosit algoritmul lui kruskal
//luam muchiile in ordine cresc a costului si le utilizam doar daca nodurile nu sunt
//in aceiasi multime (acelasi arbore partial, deoarece la un moment putem avea mai multi arbori ce vor fi reunitit)
//si reunim multimile cu tehnica de la
//disjoint
void graf::apm()
{
    unsigned int x,y;
    int cost_total=0;
    f>>this->n>>this->m;
    int tata[n+1];
    int inaltime[n+1];
    vector<pair<int,pair<int,int>>> lista_muchii_greutati;
    lista_muchii_greutati.resize(this->n+1);
    vector<pair<int,int>> sol;
    for (unsigned int i=0; i<this->m; ++i)
    {
        int cost;
        f>>x>>y>>cost;
        lista_muchii_greutati.push_back(make_pair(cost,make_pair(x,y)));
    }
    for (unsigned int i=0; i<=n; ++i)
    {
        tata[i]=0;
        inaltime[i]=1;
    }
    sort(lista_muchii_greutati.begin(),lista_muchii_greutati.end());
    unsigned int nrm=0;
    for (unsigned int i=0; i<lista_muchii_greutati.size() && nrm!=this->n-1; ++i)
    {
        int x=lista_muchii_greutati[i].second.first,y=lista_muchii_greutati[i].second.second;
        int cost_muchie=lista_muchii_greutati[i].first;
        int rx=x,ry=y,aux;
        while (tata[rx]!=0)
            rx=tata[rx];
        while (tata[ry]!=0)
            ry=tata[ry];
        if (rx==ry)
        {

        }
        else
        {

            sol.push_back(make_pair(x,y));
            ++nrm;
            cost_total+=cost_muchie;
            int hx=inaltime[rx];
            int hy=inaltime[ry];
            if (hx<hy)
            {
                tata[rx]=ry;
                inaltime[ry]=max(hy,hx+1);
            }
            else
            {
                tata[ry]=rx;
                inaltime[rx]=max(hx,hy+1);
            }
        }
        while (x!=rx)
        {
            aux=tata[x];
            tata[x]=rx;
            x=aux;
        }
        while (y!=ry)
        {
            aux=tata[y];
            tata[y]=ry;
            y=aux;
        }
    }
    g<<cost_total<<"\n"<<this->n-1<<"\n";
    for (unsigned int i=0; i<sol.size(); ++i)
        g<<sol[i].first<<" "<<sol[i].second<<"\n";
}

//se spune daca 2 elemente fac parte din aceiasi multime sau se reunesc 2 multimi
//reprezentam multimile ca arbori
//aplicam 2 imbunatatiri
//reunim arborele cu inaltime mica la cel cu inaltime mai mare (conectam radacinile)
//cand interogam tatal unui nod conectam toate nodurile vizitate direct cu radacina
//(amortizeaza cautarea tatalui la O(1))
void graf::disjoint()
{
    int n,m;
    f>>n>>m;
    int tata[n+1];
    int inaltime[n+1];
    for (unsigned int i=0; i<=n; ++i)
    {
        tata[i]=0;
        inaltime[i]=1;
    }
    int op,x,y;
    while (f>>op>>x>>y)
    {
        if (op==1)
        {
            while (tata[x]!=0)
                x=tata[x];
            while (tata[y]!=0)
                y=tata[y];
            int hx=inaltime[x];
            int hy=inaltime[y];
            if (hx<hy)
            {
                tata[x]=y;
                inaltime[y]=max(hy,hx+1);
            }
            else
            {
                tata[y]=x;
                inaltime[x]=max(hx,hy+1);
            }
        }
        else
        {
            int rx=x,ry=y,aux;
            while (tata[rx]!=0)
                rx=tata[rx];
            while (tata[ry]!=0)
                ry=tata[ry];
            if (rx==ry)
                g<<"DA\n";
            else
                g<<"NU\n";
            while (x!=rx)
            {
                aux=tata[x];
                tata[x]=rx;
                x=aux;
            }
            while (y!=ry)
            {
                aux=tata[y];
                tata[y]=ry;
                y=aux;
            }
        }
    }
}

//se calculeaza distanta minima de la un nod la celelalte
//se foloseste algoritmul lui prim(implementat cu set(sortat, un bst in c++))
//algoritmul lui prim formeaza un apm adaugand mereu la arbore muchia cu costul cel mai mic
//care nu formeaza ciclu
//la dijkstra putem verifica un nod de mai multe ori deoarece este posibil ca distanta minima sa se fi
//modificat
//de retinut ca operatiunea de decrease key(adica nodul era deja in heap si trb sa updatam distanta)
//este inlocuita de un erase
//si insert, delete se intampla doar cand nodul deja e in set (distanta!=inf)
void graf::dijkstra()
{
    unsigned int i;
    int d[this->n+1];
    set<pair<int,int>> set_noduri;
    for (i=1; i<=this->n; ++i)
    {
        d[i]=inf;
    }
    d[1]=0;
    set_noduri.insert(make_pair(0,1));
    while (!set_noduri.empty())
    {
        int nod=(*set_noduri.begin()).second;
        set_noduri.erase(set_noduri.begin());
        for (int i=0; i<this->lista_vecini[nod].size(); ++i)
        {
            int nod_dest=this->lista_vecini[nod][i].first;
            int cost_dest=this->lista_vecini[nod][i].second;
            if (d[nod]+cost_dest<d[nod_dest])
            {
                if (d[nod_dest]!=inf)
                {
                    set_noduri.erase(set_noduri.find(make_pair(d[nod_dest],nod_dest)));
                }
                d[nod_dest]=d[nod]+cost_dest;
                set_noduri.insert(make_pair(d[nod_dest],nod_dest));
            }

        }

    }
    for (unsigned int i=2; i<=this->n; ++i)
        if (d[i]!=inf)
            g<<d[i]<<" ";
        else
            g<<0<<" ";
}

//gaseste distanta minima de la un nod la muchii si in plus fata de dijkstra merge cu costuri negative
//aceiasi tehnica ca la dijkstra doar ca aici relaxam muchiile
//prelucram doar nodurile ale caror distanta a fost imbunatatita (celelalte ar fi prelucrate degeaba)
//daca un nod a fost in coada de mai mult de n ori inseamna ca sigur avem un ciclu (n muchii) negativ (altfel
//ar fi ineficient sa folosim un ciclu in drum)
void graf::Bellman_Ford()
{
    unsigned int i;
    f>>this->n>>this->m;
    int d[this->n+1];
    int nrCoada[this->n+1];
    bool inCoada[this->n+1];
    bool ciclu_negativ=false;
    queue<int> noduri;
    for (i=1; i<=this->n; ++i)
    {
        d[i]=inf;
        inCoada[i]=false;
        nrCoada[i]=0;
    }
    d[1]=0;
    noduri.push(1);
    nrCoada[1]=1;
    while (!noduri.empty() && !ciclu_negativ)
    {
        int nod_curent=noduri.front();
        noduri.pop();
        inCoada[nod_curent]=false;
        for (i=0; i<this->lista_vecini[nod_curent].size(); ++i)
        {
            int nod_dest=this->lista_vecini[nod_curent][i].first;
            int cost_dest=this->lista_vecini[nod_curent][i].second;
            if (d[nod_curent]+cost_dest<d[nod_dest])
            {
                d[nod_dest]=d[nod_curent]+cost_dest;
                if (!inCoada[nod_dest])
                {
                    if (nrCoada[nod_dest]>=this->n)
                    {
                        ciclu_negativ=true;
                    }
                    noduri.push(nod_dest);
                    inCoada[nod_dest]=true;
                    ++nrCoada[nod_dest];
                }
            }

        }
    }
    if (ciclu_negativ)
        g<<"Ciclu negativ!";
    else
        for (i=2; i<=this->n; ++i)
            g<<d[i]<<" ";
}

//incercam sa gasim o cale de la start la stop cu dfs cu proprietatea ca fiecare muchie din drum sa nu fie utilizata in totalitate
//crestem fluxul muchiilor din drum cu spatiul liber minim ramas si scadem din muchiile cu sensul opus (din muchiile din drum)
//cand nu mai gasim un drum algoritmul s-a terminat
//vom utiliza un graf rezidual care pentru fiecare muchie are si muchia cu sensul opus si capacitate 0
//acest lucru permite si operatii de undo in cazul in care am gasit un drum care nu e in fluxul maxim
void graf::maxflow()
{
    f>>this->n>>this->m;
    int capacitate[this->n+1][this->n+1];
    int flow[this->n+1][this->n+1];
    unsigned int i;
    unsigned int flux=0;
    for (i=1; i<=this->n; ++i)
        for (unsigned int j=1; j<=this->n; ++j)
            capacitate[i][j]=flow[i][j]=0;
    this->lista_vecini.resize(this->n+1);
    for (i=0; i<this->m; ++i)
    {
        int x,y,cap;
        f>>x>>y>>cap;
        this->lista_vecini[x].push_back(make_pair(y,0));
        this->lista_vecini[y].push_back(make_pair(x,0));
        capacitate[x][y]=cap;
    }
    int d[this->n+1];
    int t[this->n+1];
    do
    {
        //BFS
        for (unsigned int i=1; i<=this->n; ++i)
        {
            d[i]=-1;
            t[i]=0;
        }
        unsigned int index=0;
        vector <int> coada;
        coada.push_back(1);
        d[1]=0;
        while (index<coada.size())
        {
            int nod_curent=coada[index++];
            if (nod_curent!=this->n)
                for (unsigned int i=0; i<this->lista_vecini[nod_curent].size(); ++i)
                    if (d[this->lista_vecini[nod_curent][i].first]==-1 && capacitate[nod_curent][this->lista_vecini[nod_curent][i].first]>flow[nod_curent][this->lista_vecini[nod_curent][i].first])
                    {
                        coada.push_back(this->lista_vecini[nod_curent][i].first);
                        d[this->lista_vecini[nod_curent][i].first]=d[nod_curent]+1;
                        t[this->lista_vecini[nod_curent][i].first]=nod_curent;
                    }
        }
        //BFS end
        if (d[this->n]!=-1)
        {
            int vecin;
            for (i=0; i<this->lista_vecini[this->n].size(); ++i)
            {
                vecin=this->lista_vecini[this->n][i].first;
                if (d[vecin]!=-1 && capacitate[vecin][this->n]>flow[vecin][this->n])
                {
                    int vmin=inf;
                    int nod=this->n;
                    t[this->n]=vecin;
                    while (nod!=1)
                    {
                        vmin=min(vmin,capacitate[t[nod]][nod]-flow[t[nod]][nod]);
                        nod=t[nod];
                    }
                    if (vmin!=0)
                    {
                        nod=this->n;
                        while (nod!=1)
                        {
                            flow[t[nod]][nod]+=vmin;
                            flow[nod][t[nod]]-=vmin;
                            nod=t[nod];
                        }
                        flux+=vmin;
                    }
                }
            }
        }
    }
    while(d[this->n]!=-1);
    g<<flux;
}

//algoritmul care utilizeaza programare dinamica si formeaza o matrice unde m[i][j] este distanta de la nodul i la j
//incercam sa utilizam nodul k ca sa scurtam distanta intre nodurile i si j
//daca distanta de la i la k si k la j este mai mica ca distanta de la i la j
//atunci vom folosi distanta ce contine si nodul k
void graf::royfloyd()
{
    long long dist[101][101];
    unsigned int i,j,k;
    f>>this->n;
    for (i=1; i<=this->n; ++i)
        for (j=1; j<=this->n; ++j)
        {
            f>>dist[i][j];
            if (dist[i][j]==0)
                dist[i][j]=inf;
        }

    for (k=1; k<=this->n; ++k)
        for (i=1; i<=this->n; ++i)
            for (j=1; j<=this->n; ++j)
                if (dist[i][k]+dist[k][j]<dist[i][j])
                    dist[i][j]=dist[i][k]+dist[k][j];
    for (i=1; i<=this->n; ++i)
    {
        for (j=1; j<=this->n; ++j)
            if (dist[i][j]==inf || i==j)
                g<<0<<" ";
            else
                g<<dist[i][j]<<" ";
        g<<"\n";
    }
}

//afla diametrul unui abrore (lantul de lungime maxima intre oricare 2 noduri din arbore)
//aflam intai cel mai indepartat nod de un nod oarecare (BFS), proprietatile arborelui garanteaza ca cel mai indepartat nod e capatul lantului
//de lungime maxima
//ca sa gasim si celalalt capat facem BFS din capatul gasit anterior
void graf::darb()
{
    int d[this->n+1];
    unsigned int i;
    for (i=1; i<=this->n; ++i)
        d[i]=-1;
    BFS(1,d);
    int vmax=0;
    int nod_max=0;
    for (i=1; i<=this->n; ++i)
    {
        if (d[i]>vmax)
        {

            vmax=d[i];
            nod_max=i;
        }
        d[i]=-1;
    }
    BFS(nod_max,d);
    vmax=d[1];
    for (i=2; i<=this->n; ++i)
    {
        if (d[i]>vmax)
            vmax=d[i];
    }
    g<<vmax+1;
}

//determina daca un graf e eulerian si daca este determina un ciclu eulerian
//ne folosim de proprietatea ca un ciclu euler poate fi format din reuniunea ciclurilor disjuncte
//daca graful este eulerian:
//facem o parcurgere dfs doar ca folosim muchiile (si le stergem cand le folosim, in loc de viz ca sa nu mai folosim memorie suplimentara) deoarece este multigraf
//cand nu mai avem muchii afisam nodul
void graf::ciclueuler()
{
    f>>this->n>>this->m;
    vector<vector<int>> lista_muchii;
    unsigned int i;
    lista_muchii.resize(this->n+1);
    int n1[this->m],n2[this->m];
    bool viz[this->m];
    for (i=0; i<this->m; ++i)
    {
        int x,y;
        f>>x>>y;
        lista_muchii[x].push_back(i);
        lista_muchii[y].push_back(i);
        n1[i]=x;
        n2[i]=y;
        viz[i]=false;
    }
    for (i=1; i<=this->n; ++i)
        if (lista_muchii[i].size()%2!=0 || lista_muchii[i].size()==0)
        {
            g<<"-1";
            return;
        }
    stack <int> stiva;
    vector <int> sol;
    stiva.push(1);
    while (!stiva.empty())
    {
        int nod=stiva.top();
        if (lista_muchii[nod].size()>0)
        {
            int muchie=lista_muchii[nod].back();
            lista_muchii[nod].pop_back();
            if (!viz[muchie])
            {
                int vecin;
                if (nod!=n1[muchie])
                    vecin=n1[muchie];
                else
                    vecin=n2[muchie];
                viz[muchie]=true;
                stiva.push(vecin);
            }
        }
        else
        {
            stiva.pop();
            sol.push_back(nod);
        }
    }
    for (i=0; i<sol.size()-1; ++i)
        g<<sol[i]<<" ";
}

//determina daca un graf este hamiltonian si daca da calculeaza costul minim al unui ciclu hamiltonian
//rezolvam cu programare dinamica
//sol[i][j]=costul minim al lantului de la nodul 0 la nodul j care trece prin nodurile cu 1 din reprezentarea binara a lui i
//retinem muchiile invers pt a ne ajuta la calculul costului
//pt fiecare configuratie de noduri calculam costul minim al lantului hamiltonian pana la fiecare nod din configuratie
//pt a obtine solutia finala adaugam ultima muchie necesara pt a forma cliclul (vecinii lui 0, adica nodurile care au muchii catre 0)
void graf::cicluhamilton()
{
    f>>this->n>>this->m;
    int cost[this->n][this->n];
    int sol[1<<this->n][this->n];
    int vecin;
    unsigned int i,j;
    for (i=0; i<this->n; ++i)
        for (j=0; j<this->n; ++j)
            cost[i][j]=inf;
    for (i=0; i<(1<<this->n); ++i)
        for (j=0; j<this->n; ++j)
            sol[i][j]=inf;
    this->lista_vecini.resize(n);
    for (i=0; i<this->m; ++i)
    {
        int x,y,c;
        f>>x>>y>>c;
        this->lista_vecini[y].push_back(make_pair(x,0));
        cost[x][y]=c;
    }
    sol[1][0]=0;
    for (i=0; i<1<<this->n; ++i)
        for (j=0; j<this->n; ++j)
            if (i && (1<<j))
                for (unsigned int k=0; k<this->lista_vecini[j].size(); ++k)
                {
                    vecin=this->lista_vecini[j][k].first;
                    if (i && (1<<vecin))
                        sol[i][j]=min(sol[i][j],sol[i^(1<<j)][vecin]+cost[vecin][j]);
                }
    int cmin=inf;
    for (i=0; i<this->lista_vecini[0].size(); ++i)
    {
        vecin=this->lista_vecini[0][i].first;
        cmin=min(cmin,sol[(1<<n)-1][vecin]+cost[vecin][0]);
    }
    if (cmin!=inf)
        g<<cmin;
    else
        g<<"Nu exista solutie";

}

//bfs cuplaj construieste augumenting path
bool graf::bfs_cuplaj(int st[],int dr[],int dist[])
{
    queue<int> coada;
    unsigned int i;
    for (i=1; i<=this->n; ++i)
        if (st[i]==0)
        {
            dist[i]=0;
            coada.push(i);
        }
        else
            dist[i]=inf;
    dist[0]=inf;
    while (!coada.empty())
    {
        int nod_curent=coada.front();
        coada.pop();
        if (dist[nod_curent]<dist[0])
            for (i=0; i<this->lista_vecini[nod_curent].size(); ++i)
            {
                int vecin=this->lista_vecini[nod_curent][i].first;
                if (dist[dr[vecin]]==inf)
                {
                    dist[dr[vecin]]=dist[nod_curent]+1;
                    coada.push(dr[vecin]);
                }
            }
    }
    return dist[0]!=inf;
}

//dfs(recursiv) pt cuplaj, reconstruieste augumenting path pornind de la un nod (daca este posibil) si interschimba muchiile aflate in cuplaj cu cele care nu sunt incluse
bool graf::dfs_cuplaj(unsigned int nod,int st[], int dr[], int dist[])
{
    unsigned int i;
    if (nod!=0)
    {
        for (i=0; i<this->lista_vecini[nod].size(); ++i)
        {
            int vecin=this->lista_vecini[nod][i].first;
            if (dist[dr[vecin]]==dist[nod]+1)
                if (dfs_cuplaj(dr[vecin],st,dr,dist))
                {
                    dr[vecin]=nod;
                    st[nod]=vecin;
                    return true;
                }
        }
        dist[nod]=inf;
        return false;
    }
    return true;
}

//calculeaza cuplajul maxim intr-un graf bipartit
//afiseaza numarul de muchii dar si muchiile cuplajului
//st[i] vecinul din dreapta al nodului i(nodul i e din stanga)
//dr[i] vecinul din stanga al nodului i(nodul i e din dreapta)
//cat timp mai gasim un augumenting path mai putem creste cuplajul maxim
//augumenting path-> drum la 2 noduri care nu sunt incluse in cuplaj alternand intre muchiile folosite in cuplajul calculat la momentul actual
//si muchii libere
void graf::cuplaj()
{
    unsigned int nrm;
    unsigned int i;
    f>>this->n>>this->m>>nrm;
    this->lista_vecini.resize(this->n+1);
    for (i=0; i<nrm; ++i)
    {
        unsigned int x,y;
        f>>x>>y;
        this->lista_vecini[x].push_back(make_pair(y,0));
    }
    int st[this->n+1],dr[this->m+1],dist[this->n+1];
    for (i=0; i<=this->n+1; ++i)
        st[i]=0;
    for (i=0; i<=this->m+1; ++i)
        dr[i]=0;
    int sol=0;
    while (bfs_cuplaj(st,dr,dist))
    {
        for (i=1; i<=this->n; ++i)
            if (st[i]==0 && dfs_cuplaj(i,st,dr,dist))
                ++sol;
    }
    g<<sol<<"\n";
    for (i=0; i<=this->n; ++i)
        if (st[i])
            g<<i<<" "<<st[i]<<"\n";
}

/* format pt a adauga o problema
!schimbam #define PROBLEMA y de sus in #define PROBLEMA x!
 case x(nr unic):
    {
        f.open ("nume_fisier_intrare", std::ifstream::in);
        g.open ("nume_fisier_intrare", std::ifstream::out);
        graf nume_graf(x) (cifra optionala,vezi comentariile de la constructor);
        nume_graf.Metoda();
        break;
    }
*/

int main()
{
    switch (PROBLEMA)
    {
    case 1:
    {
        f.open ("bfs.in", std::ifstream::in);
        g.open ("bfs.out", std::ifstream::out);
        graf a(1);
        a.StartBFS();
        break;
    }
    case 2:
    {
        f.open ("dfs.in", std::ifstream::in);
        g.open ("dfs.out", std::ifstream::out);
        graf a(2);
        a.StartDFS();
        break;
    }
    case 3:
    {
        f.open("biconex.in",std::fstream::in);
        g.open("biconex.out",std::fstream::out);
        graf a(2);
        a.Startbiconex();
        break;
    }
    case 4:
    {
        f.open("ctc.in",std::fstream::in);
        g.open("ctc.out",std::fstream::out);
        graf a(3);
        a.Startctc();
        break;
    }
    case 5:
    {
        f.open("sortaret.in",std::fstream::in);
        g.open("sortaret.out",std::fstream::out);
        graf a(3);
        a.Startsortaret();
        break;
    }
    case 6:
    {
        f.open("havelhakimi.in",std::fstream::in);
        g.open("havelhakimi.out",std::fstream::out);
        graf a;
        a.havelhakimi();
        break;
    }
    case 7:
    {
        f.open("criticalConnections.in",std::fstream::in);
        g.open("criticalConnections.out",std::fstream::out);
        graf a;
        a.StartcriticalConnections();
        break;
    }
    case 8:
    {
        f.open("apm.in",std::fstream::in);
        g.open("apm.out",std::fstream::out);
        graf a;
        a.apm();
        break;
    }
    case 9:
    {
        f.open("disjoint.in",std::fstream::in);
        g.open("disjoint.out",std::fstream::out);
        graf a;
        a.disjoint();
        break;
    }
    case 10:
    {
        f.open("dijkstra.in",std::fstream::in);
        g.open("dijkstra.out",std::fstream::out);
        graf a(4);
        a.dijkstra();
        break;
    }
    case 11:
    {
        f.open("bellmanford.in",std::fstream::in);
        g.open("bellmanford.out",std::fstream::out);
        graf a(4);
        a.Bellman_Ford();
        break;
    }
    case 12:
    {
        f.open("maxflow.in",std::fstream::in);
        g.open("maxflow.out",std::fstream::out);
        graf a;
        a.maxflow();
        break;
    }
    case 13:
    {
        f.open("royfloyd.in",std::fstream::in);
        g.open("royfloyd.out",std::fstream::out);
        graf a;
        a.royfloyd();
        break;
    }
    case 14:
    {
        f.open("darb.in",std::fstream::in);
        g.open("darb.out",std::fstream::out);
        graf a(6);
        a.darb();
        break;
    }
    case 15:
    {
        f.open("ciclueuler.in",std::fstream::in);
        g.open("ciclueuler.out",std::fstream::out);
        graf a;
        a.ciclueuler();
        break;
    }
    case 16:
    {
        f.open("hamilton.in",std::fstream::in);
        g.open("hamilton.out",std::fstream::out);
        graf a;
        a.cicluhamilton();
        break;
    }
    case 17:
    {
        f.open("cuplaj.in",std::fstream::in);
        g.open("cuplaj.out",std::fstream::out);
        graf a;
        a.cuplaj();
        break;
    }
    default:
        break;
    }
    f.close();
    g.close();
}
