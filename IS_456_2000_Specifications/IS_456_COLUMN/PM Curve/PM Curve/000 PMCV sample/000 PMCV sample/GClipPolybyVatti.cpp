// GClipPolybyVatti.cpp: implementation of the GClipPolybyVatti class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#define __INC_GRENDERENGINE_H__
  #include "GRenderAll.h"
#include "GClipPolybyVatti.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//===================================================================================//
//
//                    Implement Vatti's Polygon Clipping Algolithm
//
//===================================================================================//

#define LEFT               0
#define RIGHT              1

#define ABOVE              0
#define BELOW              1

#define CLIP               0
#define SUBJ               1

#define INVERT_TRISTRIPS   FALSE


//------------------------------------------------------------------------------
//                                 MACROS
//------------------------------------------------------------------------------

#define EQ(a, b)           (fabs((a) - (b)) <= GPC_EPSILON)

#define PREV_INDEX(i, n)   ((i - 1 + n) % n)
#define NEXT_INDEX(i, n)   ((i + 1    ) % n)

#define OPTIMAL(v, i, n)   ((v[PREV_INDEX(i, n)].y != v[i].y) || \
                            (v[NEXT_INDEX(i, n)].y != v[i].y))

#define FWD_MIN(v, i, n)   ((v[PREV_INDEX(i, n)].vertex.y >= v[i].vertex.y) \
                         && (v[NEXT_INDEX(i, n)].vertex.y > v[i].vertex.y))

#define NOT_FMAX(v, i, n)   (v[NEXT_INDEX(i, n)].vertex.y > v[i].vertex.y)

#define REV_MIN(v, i, n)   ((v[PREV_INDEX(i, n)].vertex.y > v[i].vertex.y) \
                         && (v[NEXT_INDEX(i, n)].vertex.y >= v[i].vertex.y))

#define NOT_RMAX(v, i, n)   (v[PREV_INDEX(i, n)].vertex.y > v[i].vertex.y)

#define VERTEX(e,p,s,x,y)  {add_vertex(&((e)->outp[(p)]->v[(s)]), x, y); \
                            (e)->outp[(p)]->active++;}

#define P_EDGE(d,e,p,i,j)  {(d)= (e); \
                            do {(d)= (d)->prev;} while (!(d)->outp[(p)]); \
                            (i)= (d)->bot.x + (d)->dx * ((j)-(d)->bot.y);}

#define N_EDGE(d,e,p,i,j)  {(d)= (e); \
                            do {(d)= (d)->next;} while (!(d)->outp[(p)]); \
                            (i)= (d)->bot.x + (d)->dx * ((j)-(d)->bot.y);}

#define MY_NEW(p, b, s, t) {if ((b) > 0) { p= (t*) new t [b];}}

#define MY_FREE(p)         {if (p) { delete [](p); (p)= NULL;}}


/*
---------------------------------------------------------------------------
                            Private Data Types
---------------------------------------------------------------------------
*/

typedef enum                        /* Edge intersection classes         */
{
  NUL,                              /* Empty non-intersection            */
  EMX,                              /* External maximum                  */
  ELI,                              /* External left intermediate        */
  TED,                              /* Top edge                          */
  ERI,                              /* External right intermediate       */
  REDG,                             /* Right edge                        */
  IMM,                              /* Internal maximum and minimum      */
  IMN,                              /* Internal minimum                  */
  EMN,                              /* External minimum                  */
  EMM,                              /* External maximum and minimum      */
  LED,                              /* Left edge                         */
  ILI,                              /* Internal left intermediate        */
  BED,                              /* Bottom edge                       */
  IRI,                              /* Internal right intermediate       */
  IMX,                              /* Internal maximum                  */
  FUL                               /* Full non-intersection             */
} vertex_type;

typedef enum                        /* Horizontal edge states            */
{
  NH,                               /* No horizontal edge                */
  BH,                               /* Bottom horizontal edge            */
  TH                                /* Top horizontal edge               */
} h_state;

typedef enum                        /* Edge bundle state                 */
{
  UNBUNDLED,                        /* Isolated edge not within a bundle */
  BUNDLE_HEAD,                      /* Bundle head node                  */
  BUNDLE_TAIL                       /* Passive bundle tail node          */
} bundle_state;

typedef struct v_shape              /* Internal vertex list datatype     */
{
  double              x;            /* X coordinate component            */
  double              y;            /* Y coordinate component            */
  struct v_shape     *next;         /* Pointer to next vertex in list    */
} vertex_node;

typedef struct p_shape              /* Internal contour / tristrip type  */
{
  int                 active;       /* Active flag / vertex count        */
  int                 hole;         /* Hole / external contour flag      */
  vertex_node        *v[2];         /* Left and right vertex list ptrs   */
  struct p_shape     *next;         /* Pointer to next polygon contour   */
  struct p_shape     *proxy;        /* Pointer to actual structure used  */
} polygon_node;

typedef struct edge_shape
{
  gpc_vertex          vertex;       /* Piggy-backed contour vertex data  */
  gpc_vertex          bot;          /* Edge lower (x, y) coordinate      */
  gpc_vertex          top;          /* Edge upper (x, y) coordinate      */
  double              xb;           /* Scanbeam bottom x coordinate      */
  double              xt;           /* Scanbeam top x coordinate         */
  double              dx;           /* Change in x for a unit y increase */
  int                 type;         /* Clip / subject edge flag          */
  int                 bundle[2][2]; /* Bundle edge flags                 */
  int                 bside[2];     /* Bundle left / right indicators    */
  bundle_state        bstate[2];    /* Edge bundle state                 */
  polygon_node       *outp[2];      /* Output polygon / tristrip pointer */
  struct edge_shape  *prev;         /* Previous edge in the AET          */
  struct edge_shape  *next;         /* Next edge in the AET              */
  struct edge_shape  *pred;         /* Edge connected at the lower end   */
  struct edge_shape  *succ;         /* Edge connected at the upper end   */
  struct edge_shape  *next_bound;   /* Pointer to next bound in LMT      */
} edge_node;

typedef struct lmt_shape            /* Local minima table                */
{
  double              y;            /* Y coordinate at local minimum     */
  edge_node          *first_bound;  /* Pointer to bound list             */
  struct lmt_shape   *next;         /* Pointer to next local minimum     */
} lmt_node;

typedef struct sbt_t_shape          /* Scanbeam tree                     */
{
  double              y;            /* Scanbeam node y value             */
  struct sbt_t_shape *less;         /* Pointer to nodes with lower y     */
  struct sbt_t_shape *more;         /* Pointer to nodes with higher y    */
} sb_tree;

typedef struct it_shape             /* Intersection table                */
{
  edge_node          *ie[2];        /* Intersecting edge (bundle) pair   */
  gpc_vertex          point;        /* Point of intersection             */
  struct it_shape    *next;         /* The next intersection table node  */
} it_node;

typedef struct st_shape             /* Sorted edge table                 */
{
  edge_node          *edge;         /* Pointer to AET edge               */
  double              xb;           /* Scanbeam bottom x coordinate      */
  double              xt;           /* Scanbeam top x coordinate         */
  double              dx;           /* Change in x for a unit y increase */
  struct st_shape    *prev;         /* Previous edge in sorted list      */
} st_node;

typedef struct bbox_shape           /* Contour axis-aligned bounding box */
{
  double             xmin;          /* Minimum x coordinate              */
  double             ymin;          /* Minimum y coordinate              */
  double             xmax;          /* Maximum x coordinate              */
  double             ymax;          /* Maximum y coordinate              */
} bbox;


/*
===========================================================================
                               Global Data
===========================================================================
*/

/* Horizontal edge state transitions within scanbeam boundary */
const h_state next_h_state[3][6]=
{
  /*        ABOVE     BELOW     CROSS */
  /*        L   R     L   R     L   R */  
  /* NH */ {BH, TH,   TH, BH,   NH, NH},
  /* BH */ {NH, NH,   NH, NH,   TH, TH},
  /* TH */ {NH, NH,   NH, NH,   BH, BH}
};


/*
===========================================================================
                             Private Functions
===========================================================================
*/

static void reset_it(it_node **it)
{
  it_node *itn;

  while (*it)
  {
    itn= (*it)->next;
    MY_FREE(*it);
    *it= itn;
  }
}


static void reset_lmt(lmt_node **lmt)
{
  lmt_node *lmtn;

  while (*lmt)
  {
    lmtn= (*lmt)->next;
    MY_FREE(*lmt);
    *lmt= lmtn;
  }
}


static void insert_bound(edge_node **b, edge_node *e)
{
  edge_node *existing_bound;

  if (!*b)
  {
    /* Link node e to the tail of the list */
    *b= e;
  }
  else
  {
    /* Do primary sort on the x field */
    if (e[0].bot.x < (*b)[0].bot.x)
    {
      /* Insert a new node mid-list */
      existing_bound= *b;
      *b= e;
      (*b)->next_bound= existing_bound;
    }
    else
    {
      if (e[0].bot.x == (*b)[0].bot.x)
      {
        /* Do secondary sort on the dx field */
        if (e[0].dx < (*b)[0].dx)
        {
          /* Insert a new node mid-list */
          existing_bound= *b;
          *b= e;
          (*b)->next_bound= existing_bound;
        }
        else
        {
          /* Head further down the list */
          insert_bound(&((*b)->next_bound), e);
        }
      }
      else
      {
        /* Head further down the list */
        insert_bound(&((*b)->next_bound), e);
      }
    }
  }
}


static edge_node **bound_list(lmt_node **lmt, double y)
{
  lmt_node *existing_node;

  if (!*lmt)
  {
    /* Add node onto the tail end of the LMT */
    MY_NEW(*lmt, sizeof(lmt_node), "LMT insertion",lmt_node);
    (*lmt)->y= y;
    (*lmt)->first_bound= NULL;
    (*lmt)->next= NULL;
    return &((*lmt)->first_bound);
  }
  else
    if (y < (*lmt)->y)
    {
      /* Insert a new LMT node before the current node */
      existing_node= *lmt;
      MY_NEW(*lmt, sizeof(lmt_node), "LMT insertion",lmt_node);
      
      (*lmt)->y= y;
      (*lmt)->first_bound= NULL;
      (*lmt)->next= existing_node;
      return &((*lmt)->first_bound);
    }
    else
      if (y > (*lmt)->y)
        /* Head further up the LMT */
        return bound_list(&((*lmt)->next), y);
      else
        /* Use this existing LMT node */
        return &((*lmt)->first_bound);
}


static void add_to_sbtree(int *entries, sb_tree **sbtree, double y)
{
  if (!*sbtree)
  {
    /* Add a new tree node here */
    MY_NEW(*sbtree, sizeof(sb_tree), "scanbeam tree insertion",sb_tree);
    (*sbtree)->y= y;
    (*sbtree)->less= NULL;
    (*sbtree)->more= NULL;
    (*entries)++;
  }
  else
  {
    if ((*sbtree)->y > y)
    {
      /* Head into the 'less' sub-tree */
      add_to_sbtree(entries, &((*sbtree)->less), y);
    }
    else
    {
      if ((*sbtree)->y < y)
      {
        /* Head into the 'more' sub-tree */
        add_to_sbtree(entries, &((*sbtree)->more), y);
      }
    }
  }
}


static void build_sbt(int *entries, double *sbt, sb_tree *sbtree)
{
  if (sbtree->less)
    build_sbt(entries, sbt, sbtree->less);
  sbt[*entries]= sbtree->y;
  (*entries)++;
  if (sbtree->more)
    build_sbt(entries, sbt, sbtree->more);
}


static void free_sbtree(sb_tree **sbtree)
{
  if (*sbtree)
  {
    free_sbtree(&((*sbtree)->less));
    free_sbtree(&((*sbtree)->more));
    MY_FREE(*sbtree);
  }
}


static int count_optimal_vertices(gpc_vertex_list c)
{
  int result= 0;

  /* Ignore non-contributing contours */
  if (c.num_vertices > 0)
  {
    for (int i= 0; i < c.num_vertices; i++)
      /* Ignore superfluous vertices embedded in horizontal edges */
      if (OPTIMAL(c.vertex, i, c.num_vertices))
        result++;
  }
  return result;
}

static edge_node *build_lmt(lmt_node **lmt, sb_tree **sbtree,
                            int *sbt_entries, gpc_polygon *p, int type,
                            gpc_op op)
{
  int          c, min, max, num_edges, v, num_vertices;
  int          total_vertices= 0, e_index=0;
  edge_node   *e, *edge_table;

  //--------------------------------------------------------------------
  // intersect가 발생하지 않으면 gabage값이 설정되어 Crash발생!
  //
  e = NULL; edge_table = NULL;

  //--------------------------------------------------------------------

  for (c= 0; c < p->num_contours; c++)
    total_vertices+= count_optimal_vertices(p->contour[c]);

  /* Create the entire input polygon edge table in one go */
  MY_NEW(edge_table, total_vertices * sizeof(edge_node),
         "edge table creation",edge_node);

  for (c= 0; c < p->num_contours; c++)
  {
    if (p->contour[c].num_vertices < 0)
    {
      /* Ignore the non-contributing contour and repair the vertex count */
      p->contour[c].num_vertices= -p->contour[c].num_vertices;
    }
    else
    {
      /* Perform contour optimisation */
      num_vertices= 0;
      for (int i= 0; i < p->contour[c].num_vertices; i++)
        if (OPTIMAL(p->contour[c].vertex, i, p->contour[c].num_vertices))
        {
          edge_table[num_vertices].vertex.x= p->contour[c].vertex[i].x;
          edge_table[num_vertices].vertex.y= p->contour[c].vertex[i].y;

          /* Record vertex in the scanbeam table */
          add_to_sbtree(sbt_entries, sbtree,
                        edge_table[num_vertices].vertex.y);

          num_vertices++;
        }

      /* Do the contour forward pass */
      for (min= 0; min < num_vertices; min++)
      {
        /* If a forward local minimum... */
        if (FWD_MIN(edge_table, min, num_vertices))
        {
          /* Search for the next local maximum... */
          num_edges= 1;
          max= NEXT_INDEX(min, num_vertices);
          while (NOT_FMAX(edge_table, max, num_vertices))
          {
            num_edges++;
            max= NEXT_INDEX(max, num_vertices);
          }

          /* Build the next edge list */
          e= &edge_table[e_index];
          e_index+= num_edges;
          v= min;
          e[0].bstate[BELOW]= UNBUNDLED;
          e[0].bundle[BELOW][CLIP]= FALSE;
          e[0].bundle[BELOW][SUBJ]= FALSE;
          for (int i= 0; i < num_edges; i++)
          {
            e[i].xb= edge_table[v].vertex.x;
            e[i].bot.x= edge_table[v].vertex.x;
            e[i].bot.y= edge_table[v].vertex.y;

            v= NEXT_INDEX(v, num_vertices);

            e[i].top.x= edge_table[v].vertex.x;
            e[i].top.y= edge_table[v].vertex.y;
            e[i].dx= (edge_table[v].vertex.x - e[i].bot.x) /
                       (e[i].top.y - e[i].bot.y);
            e[i].type= type;
            e[i].outp[ABOVE]= NULL;
            e[i].outp[BELOW]= NULL;
            e[i].next= NULL;
            e[i].prev= NULL;
            e[i].succ= ((num_edges > 1) && (i < (num_edges - 1))) ?
                       &(e[i + 1]) : NULL;
            e[i].pred= ((num_edges > 1) && (i > 0)) ? &(e[i - 1]) : NULL;
            e[i].next_bound= NULL;
            e[i].bside[CLIP]= (op == GPC_DIFF) ? RIGHT : LEFT;
            e[i].bside[SUBJ]= LEFT;
          }
          insert_bound(bound_list(lmt, edge_table[min].vertex.y), e);
        }
      }

      /* Do the contour reverse pass */
      for (min= 0; min < num_vertices; min++)
      {
        /* If a reverse local minimum... */
        if (REV_MIN(edge_table, min, num_vertices))
        {
          /* Search for the previous local maximum... */
          num_edges= 1;
          max= PREV_INDEX(min, num_vertices);
          while (NOT_RMAX(edge_table, max, num_vertices))
          {
            num_edges++;
            max= PREV_INDEX(max, num_vertices);
          }

          /* Build the previous edge list */
          e= &edge_table[e_index];
          e_index+= num_edges;
          v= min;
          e[0].bstate[BELOW]= UNBUNDLED;
          e[0].bundle[BELOW][CLIP]= FALSE;
          e[0].bundle[BELOW][SUBJ]= FALSE;
          for (int i= 0; i < num_edges; i++)
          {
            e[i].xb= edge_table[v].vertex.x;
            e[i].bot.x= edge_table[v].vertex.x;
            e[i].bot.y= edge_table[v].vertex.y;

            v= PREV_INDEX(v, num_vertices);

            e[i].top.x= edge_table[v].vertex.x;
            e[i].top.y= edge_table[v].vertex.y;
            e[i].dx= (edge_table[v].vertex.x - e[i].bot.x) /
                       (e[i].top.y - e[i].bot.y);
            e[i].type= type;
            e[i].outp[ABOVE]= NULL;
            e[i].outp[BELOW]= NULL;
            e[i].next= NULL;
            e[i].prev= NULL;
            e[i].succ= ((num_edges > 1) && (i < (num_edges - 1))) ?
                       &(e[i + 1]) : NULL;
            e[i].pred= ((num_edges > 1) && (i > 0)) ? &(e[i - 1]) : NULL;
            e[i].next_bound= NULL;
            e[i].bside[CLIP]= (op == GPC_DIFF) ? RIGHT : LEFT;
            e[i].bside[SUBJ]= LEFT;
          }
          insert_bound(bound_list(lmt, edge_table[min].vertex.y), e);
        }
      }
    }
  }
  return edge_table;
}

static void add_edge_to_aet(edge_node **aet, edge_node *edge, edge_node *prev)
{
  if (!*aet)
  {
    /* Append edge onto the tail end of the AET */
    *aet= edge;
    edge->prev= prev;
    edge->next= NULL;
  }
  else
  {
    /* Do primary sort on the xb field */
    if (edge->xb < (*aet)->xb)
    {
      /* Insert edge here (before the AET edge) */
      edge->prev= prev;
      edge->next= *aet;
      (*aet)->prev= edge;
      *aet= edge;
    }
    else
    {
      if (edge->xb == (*aet)->xb)
      {
        /* Do secondary sort on the dx field */
        if (edge->dx < (*aet)->dx)
        {
          /* Insert edge here (before the AET edge) */
          edge->prev= prev;
          edge->next= *aet;
          (*aet)->prev= edge;
          *aet= edge;
        }
        else
        {
          /* Head further into the AET */
          add_edge_to_aet(&((*aet)->next), edge, *aet);
        }
      }
      else
      {
        /* Head further into the AET */
        add_edge_to_aet(&((*aet)->next), edge, *aet);
      }
    }
  }
}


static void add_intersection(it_node **it, edge_node *edge0, edge_node *edge1,
                             double x, double y)
{
  it_node *existing_node;

  if (!*it)
  {
    /* Append a new node to the tail of the list */
    MY_NEW(*it, sizeof(it_node), "IT insertion",it_node);
    (*it)->ie[0]= edge0;
    (*it)->ie[1]= edge1;
    (*it)->point.x= x;
    (*it)->point.y= y;
    (*it)->next= NULL;
  }
  else
  {
    if ((*it)->point.y > y)
    {
      /* Insert a new node mid-list */
      existing_node= *it;
      MY_NEW(*it, sizeof(it_node), "IT insertion",it_node);
      (*it)->ie[0]= edge0;
      (*it)->ie[1]= edge1;
      (*it)->point.x= x;
      (*it)->point.y= y;
      (*it)->next= existing_node;
    }
    else
      /* Head further down the list */
      add_intersection(&((*it)->next), edge0, edge1, x, y);
  }
}


static void add_st_edge(st_node **st, it_node **it, edge_node *edge,
                        double dy)
{
  st_node *existing_node;
  double   den, r, x, y;

  if (!*st)
  {
    /* Append edge onto the tail end of the ST */
    MY_NEW(*st, sizeof(st_node), "ST insertion",st_node);
    (*st)->edge= edge;
    (*st)->xb= edge->xb;
    (*st)->xt= edge->xt;
    (*st)->dx= edge->dx;
    (*st)->prev= NULL;
  }
  else
  {
    den= ((*st)->xt - (*st)->xb) - (edge->xt - edge->xb);

    /* If new edge and ST edge don't cross */
    if ((edge->xt >= (*st)->xt) || (edge->dx == (*st)->dx) || 
        (fabs(den) <= DBL_EPSILON))
    {
      /* No intersection - insert edge here (before the ST edge) */
      existing_node= *st;
      MY_NEW(*st, sizeof(st_node), "ST insertion",st_node);
      (*st)->edge= edge;
      (*st)->xb= edge->xb;
      (*st)->xt= edge->xt;
      (*st)->dx= edge->dx;
      (*st)->prev= existing_node;
    }
    else
    {
      /* Compute intersection between new edge and ST edge */
      r= (edge->xb - (*st)->xb) / den;
      x= (*st)->xb + r * ((*st)->xt - (*st)->xb);
      y= r * dy;

      /* Insert the edge pointers and the intersection point in the IT */
      add_intersection(it, (*st)->edge, edge, x, y);

      /* Head further into the ST */
      add_st_edge(&((*st)->prev), it, edge, dy);
    }
  }
}


static void build_intersection_table(it_node **it, edge_node *aet, double dy)
{
  st_node   *st, *stp;
  edge_node *edge;

  /* Build intersection table for the current scanbeam */
  reset_it(it);
  st= NULL;

  /* Process each AET edge */
  for (edge= aet; edge; edge= edge->next)
  {
    if ((edge->bstate[ABOVE] == BUNDLE_HEAD) ||
         edge->bundle[ABOVE][CLIP] || edge->bundle[ABOVE][SUBJ])
      add_st_edge(&st, it, edge, dy);
  }

  /* Free the sorted edge table */
  while (st)
  {
    stp= st->prev;
    MY_FREE(st);
    st= stp;
  }
}

static int count_contours(polygon_node *polygon)
{
  int          nc, nv;
  vertex_node *v, *nextv;

  for (nc= 0; polygon; polygon= polygon->next)
    if (polygon->active)
    {
      /* Count the vertices in the current contour */
      nv= 0;
      for (v= polygon->proxy->v[LEFT]; v; v= v->next)
        nv++;

      /* Record valid vertex counts in the active field */
      if (nv > 2)
      {
        polygon->active= nv;
        nc++;
      }
      else
      {
        /* Invalid contour: just free the heap */
        for (v= polygon->proxy->v[LEFT]; v; v= nextv)
        {
          nextv= v->next;
          MY_FREE(v);
        }
        polygon->active = 0;
      }
    }
  return nc;
}


static void add_left(polygon_node *p, double x, double y)
{
  vertex_node *nv;

  /* Create a new vertex node and set its fields */
  MY_NEW(nv, sizeof(vertex_node), "vertex node creation",vertex_node);
  nv->x= x;
  nv->y= y;

  /* Add vertex nv to the left end of the polygon's vertex list */
  nv->next= p->proxy->v[LEFT];

  /* Update proxy->[LEFT] to point to nv */
  p->proxy->v[LEFT]= nv;
}


static void merge_left(polygon_node *p, polygon_node *q, polygon_node *list)
{
  polygon_node *target;

  /* Label contour as a hole */
  q->proxy->hole= TRUE;

  if (p->proxy != q->proxy)
  {
    /* Assign p's vertex list to the left end of q's list */
    p->proxy->v[RIGHT]->next= q->proxy->v[LEFT];
    q->proxy->v[LEFT]= p->proxy->v[LEFT];

    /* Redirect any p->proxy references to q->proxy */
    
    for (target= p->proxy; list; list= list->next)
    {
      if (list->proxy == target)
      {
        list->active= FALSE;
        list->proxy= q->proxy;
      }
    }
  }
}


static void add_right(polygon_node *p, double x, double y)
{
  vertex_node *nv;

  /* Create a new vertex node and set its fields */
  MY_NEW(nv, sizeof(vertex_node), "vertex node creation",vertex_node);
  nv->x= x;
  nv->y= y;
  nv->next= NULL;

  /* Add vertex nv to the right end of the polygon's vertex list */
  p->proxy->v[RIGHT]->next= nv;

  /* Update proxy->v[RIGHT] to point to nv */
  p->proxy->v[RIGHT]= nv;
}


static void merge_right(polygon_node *p, polygon_node *q, polygon_node *list)
{
  polygon_node *target;

  /* Label contour as external */
  q->proxy->hole= FALSE;

  if (p->proxy != q->proxy)
  {
    /* Assign p's vertex list to the right end of q's list */
    q->proxy->v[RIGHT]->next= p->proxy->v[LEFT];
    q->proxy->v[RIGHT]= p->proxy->v[RIGHT];

    /* Redirect any p->proxy references to q->proxy */
    for (target= p->proxy; list; list= list->next)
    {
      if (list->proxy == target)
      {
        list->active= FALSE;
        list->proxy= q->proxy;
      }
    }
  }
}


static void add_local_min(polygon_node **p, edge_node *edge,
                          double x, double y)
{
  polygon_node *existing_min;
  vertex_node  *nv;

  existing_min= *p;

  MY_NEW(*p, sizeof(polygon_node), "polygon node creation",polygon_node);

  /* Create a new vertex node and set its fields */
  MY_NEW(nv, sizeof(vertex_node), "vertex node creation",vertex_node);
  nv->x= x;
  nv->y= y;
  nv->next= NULL;

  /* Initialise proxy to point to p itself */
  (*p)->proxy= (*p);
  (*p)->active= TRUE;
  (*p)->next= existing_min;

  /* Make v[LEFT] and v[RIGHT] point to new vertex nv */
  (*p)->v[LEFT]= nv;
  (*p)->v[RIGHT]= nv;

  /* Assign polygon p to the edge */
  edge->outp[ABOVE]= *p;
}


static int count_tristrips(polygon_node *tn)
{
  int total;

  for (total= 0; tn; tn= tn->next)
    if (tn->active > 2)
      total++;
  return total;
}

//------------------------------------------------------------------------------
// 
static void add_vertex(vertex_node **t, double x, double y)
{
  if (!(*t))
  {
    MY_NEW(*t, sizeof(vertex_node), "tristrip vertex creation",vertex_node);
    (*t)->x= x;
    (*t)->y= y;
    (*t)->next= NULL;
  }
  else
    /* Head further down the list */
    add_vertex(&((*t)->next), x, y);
}

static void new_tristrip(polygon_node **tn, edge_node *edge,
                         double x, double y)
{
  if (!(*tn))
  {
    MY_NEW(*tn, sizeof(polygon_node), "tristrip node creation",polygon_node);
    (*tn)->next= NULL;
    (*tn)->v[LEFT]= NULL;
    (*tn)->v[RIGHT]= NULL;
    (*tn)->active= 1;
    add_vertex(&((*tn)->v[LEFT]), x, y); 
    edge->outp[ABOVE]= *tn;
  }
  else
    /* Head further down the list */
    new_tristrip(&((*tn)->next), edge, x, y);
}

//-----------------------------------------------------------------------------
//  주어진 Polygon에 대한 Bounding Box를 생성한다. 
static bbox *create_contour_bboxes(gpc_polygon *p)
{
  bbox *box;
  int   c, v;

  MY_NEW(box, p->num_contours * sizeof(bbox), "Bounding box creation",bbox);

  /* Construct contour bounding boxes */
  for (c= 0; c < p->num_contours; c++)
  {
    /* Initialise bounding box extent */
    box[c].xmin=  DBL_MAX;
    box[c].ymin=  DBL_MAX;
    box[c].xmax= -DBL_MAX;
    box[c].ymax= -DBL_MAX;

    for (v= 0; v < p->contour[c].num_vertices; v++)
    {
      /* Adjust bounding box */
      if (p->contour[c].vertex[v].x < box[c].xmin)
        box[c].xmin= p->contour[c].vertex[v].x;
      if (p->contour[c].vertex[v].y < box[c].ymin)
        box[c].ymin= p->contour[c].vertex[v].y;
      if (p->contour[c].vertex[v].x > box[c].xmax)
        box[c].xmax= p->contour[c].vertex[v].x;
      if (p->contour[c].vertex[v].y > box[c].ymax)
          box[c].ymax= p->contour[c].vertex[v].y;
    }
  }
  return box;  
}

//-----------------------------------------------------------------------
// Min Max Test를 수행하고 현재 Operation에 부합되지 않는 Polygon 
// Vertex의 개수를 음수로 설정한다. 차후 Operation에서 제외 된다.  
static void minimax_test(gpc_polygon *subj, gpc_polygon *clip, gpc_op op)
{
  bbox *s_bbox, *c_bbox;
  int   s, c, *o_table, overlap;

  s_bbox= create_contour_bboxes(subj);
  c_bbox= create_contour_bboxes(clip);

  MY_NEW(o_table, subj->num_contours * clip->num_contours * sizeof(int),
         "overlap table creation",int);

  /* Check all subject contour bounding boxes against clip boxes */
  for (s= 0; s < subj->num_contours; s++)
    for (c= 0; c < clip->num_contours; c++)
      o_table[c * subj->num_contours + s]=
             (!((s_bbox[s].xmax < c_bbox[c].xmin) ||
                (s_bbox[s].xmin > c_bbox[c].xmax))) &&
             (!((s_bbox[s].ymax < c_bbox[c].ymin) ||
                (s_bbox[s].ymin > c_bbox[c].ymax)));

  /* For each clip contour, search for any subject contour overlaps */
  for (c= 0; c < clip->num_contours; c++)
  {
    overlap= 0;
    for (s= 0; (!overlap) && (s < subj->num_contours); s++)
      overlap= o_table[c * subj->num_contours + s];

    if (!overlap)
      /* Flag non contributing status by negating vertex count */
      clip->contour[c].num_vertices = -clip->contour[c].num_vertices;
  }  

  if (op == GPC_INT)
  {  
    /* For each subject contour, search for any clip contour overlaps */
    for (s= 0; s < subj->num_contours; s++)
    {
      overlap= 0;
      for (c= 0; (!overlap) && (c < clip->num_contours); c++)
        overlap= o_table[c * subj->num_contours + s];

      if (!overlap)
        /* Flag non contributing status by negating vertex count */
        subj->contour[s].num_vertices = -subj->contour[s].num_vertices;
    }  
  }

  MY_FREE(s_bbox);
  MY_FREE(c_bbox);
  MY_FREE(o_table);
}


/*
===========================================================================
                             Public Functions
===========================================================================
*/

void gpc_free_polygon(gpc_polygon *p)
{
  int c;

  for (c= 0; c < p->num_contours; c++)
    MY_FREE(p->contour[c].vertex);
  MY_FREE(p->hole);
  MY_FREE(p->contour);
  p->num_contours= 0;
}


void gpc_read_polygon(FILE *fp, int read_hole_flags, gpc_polygon *p)
{
  int c, v;

  fscanf_s(fp, "%d", &(p->num_contours));
  MY_NEW(p->hole, p->num_contours * sizeof(int),
         "hole flag array creation",int);
  MY_NEW(p->contour, p->num_contours
         * sizeof(gpc_vertex_list), "contour creation",gpc_vertex_list);
  for (c= 0; c < p->num_contours; c++)
  {
    fscanf_s(fp, "%d", &(p->contour[c].num_vertices));

    if (read_hole_flags)
      fscanf_s(fp, "%d", &(p->hole[c]));
    else
      p->hole[c]= FALSE; /* Assume all contours to be external */

    MY_NEW(p->contour[c].vertex, p->contour[c].num_vertices
           * sizeof(gpc_vertex), "vertex creation",gpc_vertex);
    for (v= 0; v < p->contour[c].num_vertices; v++)
      fscanf_s(fp, "%lf %lf", &(p->contour[c].vertex[v].x),
                            &(p->contour[c].vertex[v].y));
  }
}


void gpc_write_polygon(FILE *fp, int write_hole_flags, gpc_polygon *p)
{
  int c, v;

  fprintf(fp, "%d\n", p->num_contours);
  for (c= 0; c < p->num_contours; c++)
  {
    fprintf(fp, "%d\n", p->contour[c].num_vertices);

    if (write_hole_flags)
      fprintf(fp, "%d\n", p->hole[c]);
    
    for (v= 0; v < p->contour[c].num_vertices; v++)
      fprintf(fp, "% .*lf % .*lf\n",
              DBL_DIG, p->contour[c].vertex[v].x,
              DBL_DIG, p->contour[c].vertex[v].y);
  }
}


void gpc_add_contour(gpc_polygon *p, gpc_vertex_list *new_contour, int hole)
{
  int             *extended_hole, c, v;
  gpc_vertex_list *extended_contour;

  /* Create an extended hole array */
  MY_NEW(extended_hole, (p->num_contours + 1)
         * sizeof(int), "contour hole addition",int);

  /* Create an extended contour array */
  MY_NEW(extended_contour, (p->num_contours + 1)
         * sizeof(gpc_vertex_list), "contour addition",gpc_vertex_list);

  /* Copy the old contour and hole data into the extended arrays */
  for (c= 0; c < p->num_contours; c++)
  {
    extended_hole[c]= p->hole[c];
    extended_contour[c]= p->contour[c];
  }

  /* Copy the new contour and hole onto the end of the extended arrays */
  c= p->num_contours;
  extended_hole[c]= hole;
  extended_contour[c].num_vertices= new_contour->num_vertices;
  MY_NEW(extended_contour[c].vertex, new_contour->num_vertices
         * sizeof(gpc_vertex), "contour addition",gpc_vertex);
  for (v= 0; v < new_contour->num_vertices; v++)
    extended_contour[c].vertex[v]= new_contour->vertex[v];

  /* Dispose of the old contour */
  MY_FREE(p->contour);
  MY_FREE(p->hole);

  /* Update the polygon information */
  p->num_contours++;
  p->hole= extended_hole;
  p->contour= extended_contour;
}


void gpc_polygon_clip(gpc_op op, gpc_polygon *subj, gpc_polygon *clip,
                      gpc_polygon *result)
{
  sb_tree       *sbtree= NULL;
  it_node       *it= NULL, *intersect;
  edge_node     *edge, *prev_edge, *next_edge, *succ_edge, *e0, *e1;
  edge_node     *aet= NULL, *c_heap= NULL, *s_heap= NULL;
  lmt_node      *lmt= NULL, *local_min;
  polygon_node  *out_poly= NULL, *p, *q, *poly, *npoly, *cf= NULL;
  vertex_node   *vtx, *nv;
  h_state        horiz[2];
  int            in[2], exists[2], parity[2]= {LEFT, LEFT};
  int            c, v, contributing, search, scanbeam= 0, sbt_entries= 0;
  int            vclass, bl, br, tl, tr;
  double        *sbt= NULL, xb, px, yb, yt, dy, ix, iy;

  /* Test for trivial NULL result cases */
  if (((subj->num_contours == 0) && (clip->num_contours == 0))
   || ((subj->num_contours == 0) && ((op == GPC_INT) || (op == GPC_DIFF)))
   || ((clip->num_contours == 0) &&  (op == GPC_INT)))
  {
    result->num_contours= 0;
    result->hole= NULL;
    result->contour= NULL;
    return;
  }

  /* Identify potentialy contributing contours */
  if (((op == GPC_INT) || (op == GPC_DIFF))
   && (subj->num_contours > 0) && (clip->num_contours > 0))
    minimax_test(subj, clip, op);

  /* Build LMT */
  if (subj->num_contours > 0)
    s_heap= build_lmt(&lmt, &sbtree, &sbt_entries, subj, SUBJ, op);
  if (clip->num_contours > 0)
    c_heap= build_lmt(&lmt, &sbtree, &sbt_entries, clip, CLIP, op);

  /* Return a NULL result if no contours contribute */
  if (lmt == NULL)
  {
    result->num_contours= 0;
    result->hole= NULL;
    result->contour= NULL;
    reset_lmt(&lmt);
    MY_FREE(s_heap);
    MY_FREE(c_heap);
    return;
  }

  /* Build scanbeam table from scanbeam tree */
  MY_NEW(sbt, sbt_entries * sizeof(double), "sbt creation",double);
  build_sbt(&scanbeam, sbt, sbtree);
  scanbeam= 0;
  free_sbtree(&sbtree);

  /* Allow pointer re-use without causing memory leak */
  if (subj == result)
    gpc_free_polygon(subj);
  if (clip == result)
    gpc_free_polygon(clip);

  /* Invert clip polygon for difference operation */
  if (op == GPC_DIFF)
    parity[CLIP]= RIGHT;

  local_min= lmt;

  /* Process each scanbeam */
  while (scanbeam < sbt_entries)
  {
    /* Set yb and yt to the bottom and top of the scanbeam */
    yb= sbt[scanbeam++];
    if (scanbeam < sbt_entries)
    {
      yt= sbt[scanbeam];
      dy= yt - yb;
    }

    /* === SCANBEAM BOUNDARY PROCESSING ================================ */

    /* If LMT node corresponding to yb exists */
    if (local_min)
    {
      if (local_min->y == yb)
      {
        /* Add edges starting at this local minimum to the AET */
        for (edge= local_min->first_bound; edge; edge= edge->next_bound)
          add_edge_to_aet(&aet, edge, NULL);

        local_min= local_min->next;
      }
    }

    /* Set dummy previous x value */
    px= -DBL_MAX;

    /* Create bundles within AET */
    e0= aet;
    e1= aet;

    /* Set up bundle fields of first edge */
    aet->bundle[ABOVE][ aet->type]= (aet->top.y != yb);
    aet->bundle[ABOVE][!aet->type]= FALSE;
    aet->bstate[ABOVE]= UNBUNDLED;

    for (next_edge= aet->next; next_edge; next_edge= next_edge->next)
    {
      /* Set up bundle fields of next edge */
      next_edge->bundle[ABOVE][ next_edge->type]= (next_edge->top.y != yb);
      next_edge->bundle[ABOVE][!next_edge->type]= FALSE;
      next_edge->bstate[ABOVE]= UNBUNDLED;

      /* Bundle edges above the scanbeam boundary if they coincide */
      if (next_edge->bundle[ABOVE][next_edge->type])
      {
        if (EQ(e0->xb, next_edge->xb) && EQ(e0->dx, next_edge->dx)
	          && (e0->top.y != yb))
        {
          next_edge->bundle[ABOVE][ next_edge->type]^= 
            e0->bundle[ABOVE][ next_edge->type];
          next_edge->bundle[ABOVE][!next_edge->type]= 
            e0->bundle[ABOVE][!next_edge->type];
          next_edge->bstate[ABOVE]= BUNDLE_HEAD;
          e0->bundle[ABOVE][CLIP]= FALSE;
          e0->bundle[ABOVE][SUBJ]= FALSE;
          e0->bstate[ABOVE]= BUNDLE_TAIL;
        }
        e0= next_edge;
      }
    }
    
    horiz[CLIP]= NH;
    horiz[SUBJ]= NH;

    /* Process each edge at this scanbeam boundary */
    for (edge= aet; edge; edge= edge->next)
    {
      exists[CLIP]= edge->bundle[ABOVE][CLIP] + 
                   (edge->bundle[BELOW][CLIP] << 1);
      exists[SUBJ]= edge->bundle[ABOVE][SUBJ] + 
                   (edge->bundle[BELOW][SUBJ] << 1);

      if (exists[CLIP] || exists[SUBJ])
      {
        /* Set bundle side */
        edge->bside[CLIP]= parity[CLIP];
        edge->bside[SUBJ]= parity[SUBJ];

        /* Determine contributing status and quadrant occupancies */
        switch (op)
        {
        case GPC_DIFF:
        case GPC_INT:
          contributing= (exists[CLIP] && (parity[SUBJ] || horiz[SUBJ]))
                     || (exists[SUBJ] && (parity[CLIP] || horiz[CLIP]))
                     || (exists[CLIP] && exists[SUBJ]
                     && (parity[CLIP] == parity[SUBJ]));
          br= (parity[CLIP])
           && (parity[SUBJ]);
          bl= (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
           && (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
          tr= (parity[CLIP] ^ (horiz[CLIP]!=NH))
           && (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
          tl= (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP]) 
           && (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
          break;
        case GPC_XOR:
          contributing= exists[CLIP] || exists[SUBJ];
          br= (parity[CLIP])
            ^ (parity[SUBJ]);
          bl= (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
            ^ (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
          tr= (parity[CLIP] ^ (horiz[CLIP]!=NH))
            ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
          tl= (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP]) 
            ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
          break;
        case GPC_UNION:
          contributing= (exists[CLIP] && (!parity[SUBJ] || horiz[SUBJ]))
                     || (exists[SUBJ] && (!parity[CLIP] || horiz[CLIP]))
                     || (exists[CLIP] && exists[SUBJ]
                     && (parity[CLIP] == parity[SUBJ]));
          br= (parity[CLIP])
           || (parity[SUBJ]);
          bl= (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
           || (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
          tr= (parity[CLIP] ^ (horiz[CLIP]!=NH))
           || (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
          tl= (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP]) 
           || (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
          break;
        }

        /* Update parity */
        parity[CLIP]^= edge->bundle[ABOVE][CLIP];
        parity[SUBJ]^= edge->bundle[ABOVE][SUBJ];

        /* Update horizontal state */
        if (exists[CLIP])         
          horiz[CLIP]=
            next_h_state[horiz[CLIP]]
                        [((exists[CLIP] - 1) << 1) + parity[CLIP]];
        if (exists[SUBJ])         
          horiz[SUBJ]=
            next_h_state[horiz[SUBJ]]
                        [((exists[SUBJ] - 1) << 1) + parity[SUBJ]];

        vclass= tr + (tl << 1) + (br << 2) + (bl << 3);

        if (contributing)
        {
          xb= edge->xb;

          switch (vclass)
          {
          case EMN:
          case IMN:
            add_local_min(&out_poly, edge, xb, yb);
            px= xb;
            cf= edge->outp[ABOVE];
            break;
          case ERI:
            if (xb != px)
            {
              add_right(cf, xb, yb);
              px= xb;
            }
            edge->outp[ABOVE]= cf;
            cf= NULL;
            break;
          case ELI:
            add_left(edge->outp[BELOW], xb, yb);
            px= xb;
            cf= edge->outp[BELOW];
            break;
          case EMX:
            if (xb != px)
            {
              add_left(cf, xb, yb);
              px= xb;
            }
            merge_right(cf, edge->outp[BELOW], out_poly);
            cf= NULL;
            break;
          case ILI:
            if (xb != px)
            {
              add_left(cf, xb, yb);
              px= xb;
            }
            edge->outp[ABOVE]= cf;
            cf= NULL;
            break;
          case IRI:
            add_right(edge->outp[BELOW], xb, yb);
            px= xb;
            cf= edge->outp[BELOW];
            edge->outp[BELOW]= NULL;
            break;
          case IMX:
            if (xb != px)
            {
              add_right(cf, xb, yb);
              px= xb;
            }
            merge_left(cf, edge->outp[BELOW], out_poly);
            cf= NULL;
            edge->outp[BELOW]= NULL;
            break;
          case IMM:
            if (xb != px)
            {
              add_right(cf, xb, yb);
              px= xb;
            } 
            merge_left(cf, edge->outp[BELOW], out_poly);
            edge->outp[BELOW]= NULL;
            add_local_min(&out_poly, edge, xb, yb);
            cf= edge->outp[ABOVE];
            break;
          case EMM:
            if (xb != px)
            {
              add_left(cf, xb, yb);
              px= xb;
            }
            merge_right(cf, edge->outp[BELOW], out_poly);
            edge->outp[BELOW]= NULL;
            add_local_min(&out_poly, edge, xb, yb);
            cf= edge->outp[ABOVE];
            break;
          case LED:
            if (edge->bot.y == yb)
              add_left(edge->outp[BELOW], xb, yb);
            edge->outp[ABOVE]= edge->outp[BELOW];
            px= xb;
            break;
          case REDG:
            if (edge->bot.y == yb)
              add_right(edge->outp[BELOW], xb, yb);
            edge->outp[ABOVE]= edge->outp[BELOW];
            px= xb;
            break;
          default:
            break;
          } /* End of switch */
        } /* End of contributing conditional */
      } /* End of edge exists conditional */
    } /* End of AET loop */

    /* Delete terminating edges from the AET, otherwise compute xt */
    for (edge= aet; edge; edge= edge->next)
    {
      if (edge->top.y == yb)
      {
        prev_edge= edge->prev;
        next_edge= edge->next;
        if (prev_edge)
          prev_edge->next= next_edge;
        else
          aet= next_edge;
        if (next_edge)
          next_edge->prev= prev_edge;

        /* Copy bundle head state to the adjacent tail edge if required */
        if ((edge->bstate[BELOW] == BUNDLE_HEAD) && prev_edge)
        {
          if (prev_edge->bstate[BELOW] == BUNDLE_TAIL)
          {
            prev_edge->outp[BELOW]= edge->outp[BELOW];
            prev_edge->bstate[BELOW]= UNBUNDLED;
            if (prev_edge->prev)
              if (prev_edge->prev->bstate[BELOW] == BUNDLE_TAIL)
                prev_edge->bstate[BELOW]= BUNDLE_HEAD;
          }
        }
      }
      else
      {
        if (edge->top.y == yt)
          edge->xt= edge->top.x;
        else
          edge->xt= edge->bot.x + edge->dx * (yt - edge->bot.y);
      }
    }

    if (scanbeam < sbt_entries)
    {
      /* === SCANBEAM INTERIOR PROCESSING ============================== */

      build_intersection_table(&it, aet, dy);

      /* Process each node in the intersection table */
      for (intersect= it; intersect; intersect= intersect->next)
      {
        e0= intersect->ie[0];
        e1= intersect->ie[1];

        /* Only generate output for contributing intersections */
        if ((e0->bundle[ABOVE][CLIP] || e0->bundle[ABOVE][SUBJ])
         && (e1->bundle[ABOVE][CLIP] || e1->bundle[ABOVE][SUBJ]))
        {
          p= e0->outp[ABOVE];
          q= e1->outp[ABOVE];
          ix= intersect->point.x;
          iy= intersect->point.y + yb;
 
          in[CLIP]= ( e0->bundle[ABOVE][CLIP] && !e0->bside[CLIP])
                 || ( e1->bundle[ABOVE][CLIP] &&  e1->bside[CLIP])
                 || (!e0->bundle[ABOVE][CLIP] && !e1->bundle[ABOVE][CLIP]
                     && e0->bside[CLIP] && e1->bside[CLIP]);
          in[SUBJ]= ( e0->bundle[ABOVE][SUBJ] && !e0->bside[SUBJ])
                 || ( e1->bundle[ABOVE][SUBJ] &&  e1->bside[SUBJ])
                 || (!e0->bundle[ABOVE][SUBJ] && !e1->bundle[ABOVE][SUBJ]
                     && e0->bside[SUBJ] && e1->bside[SUBJ]);
       
          /* Determine quadrant occupancies */
          switch (op)
          {
          case GPC_DIFF:
          case GPC_INT:
            tr= (in[CLIP]) && (in[SUBJ]);
            tl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP]) && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
            br= (in[CLIP] ^ e0->bundle[ABOVE][CLIP]) && (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            bl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
             && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            break;
          case GPC_XOR:
            tr= (in[CLIP])^ (in[SUBJ]);
            tl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP])^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
            br= (in[CLIP] ^ e0->bundle[ABOVE][CLIP])^ (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            bl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
              ^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            break;
          case GPC_UNION:
            tr= (in[CLIP]) || (in[SUBJ]);
            tl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP]) || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
            br= (in[CLIP] ^ e0->bundle[ABOVE][CLIP]) || (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            bl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
             || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            break;
          }
	  
          vclass= tr + (tl << 1) + (br << 2) + (bl << 3);

          switch (vclass)
          {
          case EMN:
            add_local_min(&out_poly, e0, ix, iy);
            e1->outp[ABOVE]= e0->outp[ABOVE];
            break;
          case ERI:
            if (p)
            {
              add_right(p, ix, iy);
              e1->outp[ABOVE]= p;
              e0->outp[ABOVE]= NULL;
            }
            break;
          case ELI:
            if (q)
            {
              add_left(q, ix, iy);
              e0->outp[ABOVE]= q;
              e1->outp[ABOVE]= NULL;
            }
            break;
          case EMX:
            if (p && q)
            {
              add_left(p, ix, iy);
              merge_right(p, q, out_poly);
              e0->outp[ABOVE]= NULL;
              e1->outp[ABOVE]= NULL;
            }
            break;
          case IMN:
            add_local_min(&out_poly, e0, ix, iy);
            e1->outp[ABOVE]= e0->outp[ABOVE];
            break;
          case ILI:
            if (p)
            {
              add_left(p, ix, iy);
              e1->outp[ABOVE]= p;
              e0->outp[ABOVE]= NULL;
            }
            break;
          case IRI:
            if (q)
            {
              add_right(q, ix, iy);
              e0->outp[ABOVE]= q;
              e1->outp[ABOVE]= NULL;
            }
            break;
          case IMX:
            if (p && q)
            {
              add_right(p, ix, iy);
              merge_left(p, q, out_poly);
              e0->outp[ABOVE]= NULL;
              e1->outp[ABOVE]= NULL;
            }
            break;
          case IMM:
            if (p && q)
            {
              add_right(p, ix, iy);
              merge_left(p, q, out_poly);
              add_local_min(&out_poly, e0, ix, iy);
              e1->outp[ABOVE]= e0->outp[ABOVE];
            }
            break;
          case EMM:
            if (p && q)
            {
              add_left(p, ix, iy);
              merge_right(p, q, out_poly);
              add_local_min(&out_poly, e0, ix, iy);
              e1->outp[ABOVE]= e0->outp[ABOVE];
            }
            break;
          default:
            break;
          } /* End of switch */
	      } /* End of contributing intersection conditional */

        /* Swap bundle sides in response to edge crossing */
        if (e0->bundle[ABOVE][CLIP])
	        e1->bside[CLIP]= !e1->bside[CLIP];
        if (e1->bundle[ABOVE][CLIP])
	        e0->bside[CLIP]= !e0->bside[CLIP];
        if (e0->bundle[ABOVE][SUBJ])
	        e1->bside[SUBJ]= !e1->bside[SUBJ];
        if (e1->bundle[ABOVE][SUBJ])
	        e0->bside[SUBJ]= !e0->bside[SUBJ];

        /* Swap e0 and e1 bundles in the AET */
        prev_edge= e0->prev;
        next_edge= e1->next;
        if (next_edge)
          next_edge->prev= e0;

        if (e0->bstate[ABOVE] == BUNDLE_HEAD)
        {
          search= TRUE;
          while (search)
          {
            prev_edge= prev_edge->prev;
            if (prev_edge)
            {
              if (prev_edge->bstate[ABOVE] != BUNDLE_TAIL)
                search= FALSE;
            }
            else
              search= FALSE;
          }
        }
        if (!prev_edge)
        {
          aet->prev= e1;
          e1->next= aet;
          aet= e0->next;
        }
        else
        {
          prev_edge->next->prev= e1;
          e1->next= prev_edge->next;
          prev_edge->next= e0->next;
        }
        e0->next->prev= prev_edge;
        e1->next->prev= e1;
        e0->next= next_edge;
      } /* End of IT loop*/

      /* Prepare for next scanbeam */
      for (edge= aet; edge; edge= next_edge)
      {
        next_edge= edge->next;
        succ_edge= edge->succ;

        if ((edge->top.y == yt) && succ_edge)
        {
          /* Replace AET edge by its successor */
          succ_edge->outp[BELOW]= edge->outp[ABOVE];
          succ_edge->bstate[BELOW]= edge->bstate[ABOVE];
          succ_edge->bundle[BELOW][CLIP]= edge->bundle[ABOVE][CLIP];
          succ_edge->bundle[BELOW][SUBJ]= edge->bundle[ABOVE][SUBJ];
          prev_edge= edge->prev;
          if (prev_edge)
            prev_edge->next= succ_edge;
          else
            aet= succ_edge;
          if (next_edge)
            next_edge->prev= succ_edge;
          succ_edge->prev= prev_edge;
          succ_edge->next= next_edge;
        }
        else
        {
          /* Update this edge */
          edge->outp[BELOW]        = edge->outp  [ABOVE]      ;
          edge->bstate[BELOW]      = edge->bstate[ABOVE]      ;
          edge->bundle[BELOW][CLIP]= edge->bundle[ABOVE][CLIP];
          edge->bundle[BELOW][SUBJ]= edge->bundle[ABOVE][SUBJ];
          edge->xb= edge->xt;
        }
        edge->outp[ABOVE]= NULL;
      }
    }
  } /* === END OF SCANBEAM PROCESSING ================================== */

  /* Generate result polygon from out_poly */
  result->contour= NULL;
  result->hole= NULL;
  result->num_contours= count_contours(out_poly);
  if (result->num_contours > 0)
  {
    MY_NEW(result->hole, result->num_contours
           * sizeof(int), "hole flag table creation",int);
    MY_NEW(result->contour, result->num_contours
           * sizeof(gpc_vertex_list), "contour creation",gpc_vertex_list);

    c= 0;
    for (poly= out_poly; poly; poly= npoly)
    {
      npoly= poly->next;
      if (poly->active)
      {
        result->hole[c]= poly->proxy->hole;
        result->contour[c].num_vertices= poly->active;
        MY_NEW(result->contour[c].vertex,
          result->contour[c].num_vertices * sizeof(gpc_vertex),
          "vertex creation",gpc_vertex);
      
        v= result->contour[c].num_vertices - 1;
        for (vtx= poly->proxy->v[LEFT]; vtx; vtx= nv)
        {
          nv= vtx->next;
          result->contour[c].vertex[v].x= vtx->x;
          result->contour[c].vertex[v].y= vtx->y;
          MY_FREE(vtx);
          v--;
        }
        c++;
      }
      MY_FREE(poly);
    }
  }

  /* Tidy up */
  reset_it(&it);
  reset_lmt(&lmt);
  MY_FREE(c_heap);
  MY_FREE(s_heap);
  MY_FREE(sbt);
}


void gpc_free_tristrip(gpc_tristrip *t)
{
  int s;

  for (s= 0; s < t->num_strips; s++)
    MY_FREE(t->strip[s].vertex);
  MY_FREE(t->strip);
  t->num_strips= 0;
}


void gpc_polygon_to_tristrip(gpc_polygon *s, gpc_tristrip *t)
{
  gpc_polygon c;

  c.num_contours= 0;
  c.hole= NULL;
  c.contour= NULL;
  gpc_tristrip_clip(GPC_DIFF, s, &c, t);
}


void gpc_tristrip_clip(gpc_op op, gpc_polygon *subj, gpc_polygon *clip,
                       gpc_tristrip *result)
{
  sb_tree       *sbtree= NULL;
  it_node       *it= NULL, *intersect;
  edge_node     *edge, *prev_edge, *next_edge, *succ_edge, *e0, *e1;
  edge_node     *aet= NULL, *c_heap= NULL, *s_heap= NULL, *cf;
  lmt_node      *lmt= NULL, *local_min;
  polygon_node  *tlist= NULL, *tn, *tnn, *p, *q;
  vertex_node   *lt, *ltn, *rt, *rtn;
  h_state        horiz[2];
  vertex_type    cft;
  int            in[2], exists[2], parity[2]= {LEFT, LEFT};
  int            s, v, contributing, search, scanbeam= 0, sbt_entries= 0;
  int            vclass, bl, br, tl, tr;
  double        *sbt= NULL, xb, px, nx, yb, yt, dy, ix, iy;

  /* Test for trivial NULL result cases */
  if (((subj->num_contours == 0) && (clip->num_contours == 0))
   || ((subj->num_contours == 0) && ((op == GPC_INT) || (op == GPC_DIFF)))
   || ((clip->num_contours == 0) &&  (op == GPC_INT)))
  {
    result->num_strips= 0;
    result->strip= NULL;
    return;
  }

  /* Identify potentialy contributing contours */
  if (((op == GPC_INT) || (op == GPC_DIFF))
   && (subj->num_contours > 0) && (clip->num_contours > 0))
    minimax_test(subj, clip, op);

  /* Build LMT */
  if (subj->num_contours > 0)
    s_heap= build_lmt(&lmt, &sbtree, &sbt_entries, subj, SUBJ, op);
  if (clip->num_contours > 0)
    c_heap= build_lmt(&lmt, &sbtree, &sbt_entries, clip, CLIP, op);

  /* Return a NULL result if no contours contribute */
  if (lmt == NULL)
  {
    result->num_strips= 0;
    result->strip= NULL;
    reset_lmt(&lmt);
    MY_FREE(s_heap);
    MY_FREE(c_heap);
    return;
  }

  /* Build scanbeam table from scanbeam tree */
  MY_NEW(sbt, sbt_entries * sizeof(double), "sbt creation",double);
  build_sbt(&scanbeam, sbt, sbtree);
  scanbeam= 0;
  free_sbtree(&sbtree);

  /* Invert clip polygon for difference operation */
  if (op == GPC_DIFF)
    parity[CLIP]= RIGHT;

  local_min= lmt;

  /* Process each scanbeam */
  while (scanbeam < sbt_entries)
  {
    /* Set yb and yt to the bottom and top of the scanbeam */
    yb= sbt[scanbeam++];
    if (scanbeam < sbt_entries)
    {
      yt= sbt[scanbeam];
      dy= yt - yb;
    }

    /* === SCANBEAM BOUNDARY PROCESSING ================================ */

    /* If LMT node corresponding to yb exists */
    if (local_min)
    {
      if (local_min->y == yb)
      {
        /* Add edges starting at this local minimum to the AET */
        for (edge= local_min->first_bound; edge; edge= edge->next_bound)
          add_edge_to_aet(&aet, edge, NULL);

        local_min= local_min->next;
      }
    }

    /* Set dummy previous x value */
    px= -DBL_MAX;

    /* Create bundles within AET */
    e0= aet;
    e1= aet;

    /* Set up bundle fields of first edge */
    aet->bundle[ABOVE][ aet->type]= (aet->top.y != yb);
    aet->bundle[ABOVE][!aet->type]= FALSE;
    aet->bstate[ABOVE]= UNBUNDLED;

    for (next_edge= aet->next; next_edge; next_edge= next_edge->next)
    {
      /* Set up bundle fields of next edge */
      next_edge->bundle[ABOVE][ next_edge->type]= (next_edge->top.y != yb);
      next_edge->bundle[ABOVE][!next_edge->type]= FALSE;
      next_edge->bstate[ABOVE]= UNBUNDLED;

      /* Bundle edges above the scanbeam boundary if they coincide */
      if (next_edge->bundle[ABOVE][next_edge->type])
      {
        if (EQ (e0->xb, next_edge->xb) && EQ(e0->dx, next_edge->dx)
	          && (e0->top.y != yb))
        {
          next_edge->bundle[ABOVE][ next_edge->type]^= 
            e0->bundle[ABOVE][ next_edge->type];
          next_edge->bundle[ABOVE][!next_edge->type]= 
            e0->bundle[ABOVE][!next_edge->type]; 
          next_edge->bstate[ABOVE]= BUNDLE_HEAD;
          e0->bundle[ABOVE][CLIP]= FALSE;
          e0->bundle[ABOVE][SUBJ]= FALSE;
          e0->bstate[ABOVE]= BUNDLE_TAIL;
        }
        e0= next_edge;
      }
    }

    horiz[CLIP]= NH;
    horiz[SUBJ]= NH;

    /* Process each edge at this scanbeam boundary */
    for (edge= aet; edge; edge= edge->next)
    {
      exists[CLIP]= edge->bundle[ABOVE][CLIP] + 
                   (edge->bundle[BELOW][CLIP] << 1);
      exists[SUBJ]= edge->bundle[ABOVE][SUBJ] + 
                   (edge->bundle[BELOW][SUBJ] << 1);

      if (exists[CLIP] || exists[SUBJ])
      {
        /* Set bundle side */
        edge->bside[CLIP]= parity[CLIP];
        edge->bside[SUBJ]= parity[SUBJ];

        /* Determine contributing status and quadrant occupancies */
        switch (op)
        {
        case GPC_DIFF:
        case GPC_INT:
          contributing= (exists[CLIP] && (parity[SUBJ] || horiz[SUBJ]))
                     || (exists[SUBJ] && (parity[CLIP] || horiz[CLIP]))
                     || (exists[CLIP] && exists[SUBJ]
                     && (parity[CLIP] == parity[SUBJ]));
          br= (parity[CLIP])
           && (parity[SUBJ]);
          bl= (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
           && (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
          tr= (parity[CLIP] ^ (horiz[CLIP]!=NH))
           && (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
          tl= (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP]) 
           && (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
          break;
        case GPC_XOR:
          contributing= exists[CLIP] || exists[SUBJ];
          br= (parity[CLIP])
            ^ (parity[SUBJ]);
          bl= (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
            ^ (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
          tr= (parity[CLIP] ^ (horiz[CLIP]!=NH))
            ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
          tl= (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
            ^ (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
          break;
        case GPC_UNION:
          contributing= (exists[CLIP] && (!parity[SUBJ] || horiz[SUBJ]))
                     || (exists[SUBJ] && (!parity[CLIP] || horiz[CLIP]))
                     || (exists[CLIP] && exists[SUBJ]
                     && (parity[CLIP] == parity[SUBJ]));
          br= (parity[CLIP])
           || (parity[SUBJ]);
          bl= (parity[CLIP] ^ edge->bundle[ABOVE][CLIP])
           || (parity[SUBJ] ^ edge->bundle[ABOVE][SUBJ]);
          tr= (parity[CLIP] ^ (horiz[CLIP]!=NH))
           || (parity[SUBJ] ^ (horiz[SUBJ]!=NH));
          tl= (parity[CLIP] ^ (horiz[CLIP]!=NH) ^ edge->bundle[BELOW][CLIP])
           || (parity[SUBJ] ^ (horiz[SUBJ]!=NH) ^ edge->bundle[BELOW][SUBJ]);
          break;
        }

        /* Update parity */
        parity[CLIP]^= edge->bundle[ABOVE][CLIP];
        parity[SUBJ]^= edge->bundle[ABOVE][SUBJ];

        /* Update horizontal state */
        if (exists[CLIP])         
          horiz[CLIP]=
            next_h_state[horiz[CLIP]]
                        [((exists[CLIP] - 1) << 1) + parity[CLIP]];
        if (exists[SUBJ])         
          horiz[SUBJ]=
            next_h_state[horiz[SUBJ]]
                        [((exists[SUBJ] - 1) << 1) + parity[SUBJ]];
        
        vclass= tr + (tl << 1) + (br << 2) + (bl << 3);

        if (contributing)
        {
          xb= edge->xb;

          switch (vclass)
          {
          case EMN:
            new_tristrip(&tlist, edge, xb, yb);
            cf= edge;
            break;
          case ERI:
            edge->outp[ABOVE]= cf->outp[ABOVE];
            if (xb != cf->xb)
              VERTEX(edge, ABOVE, RIGHT, xb, yb);
            cf= NULL;
            break;
          case ELI:
            VERTEX(edge, BELOW, LEFT, xb, yb);
            edge->outp[ABOVE]= NULL;
            cf= edge;
            break;
          case EMX:
            if (xb != cf->xb)
              VERTEX(edge, BELOW, RIGHT, xb, yb);
            edge->outp[ABOVE]= NULL;
            cf= NULL;
            break;
          case IMN:
            if (cft == LED)
            {
              if (cf->bot.y != yb)
                VERTEX(cf, BELOW, LEFT, cf->xb, yb);
              new_tristrip(&tlist, cf, cf->xb, yb);
            }
            edge->outp[ABOVE]= cf->outp[ABOVE];
            VERTEX(edge, ABOVE, RIGHT, xb, yb);
            break;
          case ILI:
            new_tristrip(&tlist, edge, xb, yb);
            cf= edge;
            cft= ILI;
            break;
          case IRI:
            if (cft == LED)
            {
              if (cf->bot.y != yb)
                VERTEX(cf, BELOW, LEFT, cf->xb, yb);
              new_tristrip(&tlist, cf, cf->xb, yb);
            }
            VERTEX(edge, BELOW, RIGHT, xb, yb);
            edge->outp[ABOVE]= NULL;
            break;
          case IMX:
            VERTEX(edge, BELOW, LEFT, xb, yb);
            edge->outp[ABOVE]= NULL;
            cft= IMX;
            break;
	        case IMM:
            VERTEX(edge, BELOW, LEFT, xb, yb);
            edge->outp[ABOVE]= cf->outp[ABOVE];
            if (xb != cf->xb)
              VERTEX(cf, ABOVE, RIGHT, xb, yb);
            cf= edge;
            break;
          case EMM:
            VERTEX(edge, BELOW, RIGHT, xb, yb);
            edge->outp[ABOVE]= NULL;
            new_tristrip(&tlist, edge, xb, yb);
            cf= edge;
            break;
          case LED:
            if (edge->bot.y == yb)
              VERTEX(edge, BELOW, LEFT, xb, yb);
            edge->outp[ABOVE]= edge->outp[BELOW];
            cf= edge;
            cft= LED;
            break;
          case REDG:
            edge->outp[ABOVE]= cf->outp[ABOVE];
            if (cft == LED)
            {
              if (cf->bot.y == yb)
              {
                VERTEX(edge, BELOW, RIGHT, xb, yb);
              }
              else
              {
                if (edge->bot.y == yb)
                {
                  VERTEX(cf, BELOW, LEFT, cf->xb, yb);
                  VERTEX(edge, BELOW, RIGHT, xb, yb);
                }
              }
            }
            else
            {
              VERTEX(edge, BELOW, RIGHT, xb, yb);
              VERTEX(edge, ABOVE, RIGHT, xb, yb);
            }
            cf= NULL;
            break;
          default:
            break;
          } /* End of switch */
        } /* End of contributing conditional */
      } /* End of edge exists conditional */
    } /* End of AET loop */

    /* Delete terminating edges from the AET, otherwise compute xt */
    for (edge= aet; edge; edge= edge->next)
    {
      if (edge->top.y == yb)
      {
        prev_edge= edge->prev;
        next_edge= edge->next;
        if (prev_edge)
          prev_edge->next= next_edge;
        else
          aet= next_edge;
        if (next_edge)
          next_edge->prev= prev_edge;

        /* Copy bundle head state to the adjacent tail edge if required */
        if ((edge->bstate[BELOW] == BUNDLE_HEAD) && prev_edge)
        {
          if (prev_edge->bstate[BELOW] == BUNDLE_TAIL)
          {
            prev_edge->outp[BELOW]= edge->outp[BELOW];
            prev_edge->bstate[BELOW]= UNBUNDLED;
            if (prev_edge->prev)
              if (prev_edge->prev->bstate[BELOW] == BUNDLE_TAIL)
                prev_edge->bstate[BELOW]= BUNDLE_HEAD;
          }
        }
      }
      else
      {
        if (edge->top.y == yt)
          edge->xt= edge->top.x;
        else
          edge->xt= edge->bot.x + edge->dx * (yt - edge->bot.y);
      }
    }

    if (scanbeam < sbt_entries)
    {
      /* === SCANBEAM INTERIOR PROCESSING ============================== */
  
      build_intersection_table(&it, aet, dy);

      /* Process each node in the intersection table */
      for (intersect= it; intersect; intersect= intersect->next)
      {
        e0= intersect->ie[0];
        e1= intersect->ie[1];

        /* Only generate output for contributing intersections */
        if ((e0->bundle[ABOVE][CLIP] || e0->bundle[ABOVE][SUBJ])
         && (e1->bundle[ABOVE][CLIP] || e1->bundle[ABOVE][SUBJ]))
        {
          p= e0->outp[ABOVE];
          q= e1->outp[ABOVE];
          ix= intersect->point.x;
          iy= intersect->point.y + yb;

          in[CLIP]= ( e0->bundle[ABOVE][CLIP] && !e0->bside[CLIP])
                 || ( e1->bundle[ABOVE][CLIP] &&  e1->bside[CLIP])
                 || (!e0->bundle[ABOVE][CLIP] && !e1->bundle[ABOVE][CLIP]
                     && e0->bside[CLIP] && e1->bside[CLIP]);
          in[SUBJ]= ( e0->bundle[ABOVE][SUBJ] && !e0->bside[SUBJ])
                 || ( e1->bundle[ABOVE][SUBJ] &&  e1->bside[SUBJ])
                 || (!e0->bundle[ABOVE][SUBJ] && !e1->bundle[ABOVE][SUBJ]
                     && e0->bside[SUBJ] && e1->bside[SUBJ]);

          /* Determine quadrant occupancies */
          switch (op)
          {
          case GPC_DIFF:
          case GPC_INT:
            tr= (in[CLIP])
             && (in[SUBJ]);
            tl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
             && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
            br= (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
             && (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            bl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
             && (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            break;
          case GPC_XOR:
            tr= (in[CLIP])
              ^ (in[SUBJ]);
            tl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
              ^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
            br= (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
              ^ (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            bl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
              ^ (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            break;
          case GPC_UNION:
            tr= (in[CLIP])
             || (in[SUBJ]);
            tl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP])
             || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ]);
            br= (in[CLIP] ^ e0->bundle[ABOVE][CLIP])
             || (in[SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            bl= (in[CLIP] ^ e1->bundle[ABOVE][CLIP] ^ e0->bundle[ABOVE][CLIP])
             || (in[SUBJ] ^ e1->bundle[ABOVE][SUBJ] ^ e0->bundle[ABOVE][SUBJ]);
            break;
          }

          vclass= tr + (tl << 1) + (br << 2) + (bl << 3);

          switch (vclass)
          {
          case EMN:
            new_tristrip(&tlist, e1, ix, iy);
            e0->outp[ABOVE]= e1->outp[ABOVE];
            break;
          case ERI:
            if (p)
            {
              P_EDGE(prev_edge, e0, ABOVE, px, iy);
              VERTEX(prev_edge, ABOVE, LEFT, px, iy);
              VERTEX(e0, ABOVE, RIGHT, ix, iy);
              e1->outp[ABOVE]= e0->outp[ABOVE];
              e0->outp[ABOVE]= NULL;
            }
            break;
          case ELI:
            if (q)
            {
              N_EDGE(next_edge, e1, ABOVE, nx, iy);
              VERTEX(e1, ABOVE, LEFT, ix, iy);
              VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
              e0->outp[ABOVE]= e1->outp[ABOVE];
              e1->outp[ABOVE]= NULL;
            }
            break;
          case EMX:
            if (p && q)
            {
              VERTEX(e0, ABOVE, LEFT, ix, iy);
              e0->outp[ABOVE]= NULL;
              e1->outp[ABOVE]= NULL;
            }
            break;
          case IMN:
            P_EDGE(prev_edge, e0, ABOVE, px, iy);
            VERTEX(prev_edge, ABOVE, LEFT, px, iy);
            N_EDGE(next_edge, e1, ABOVE, nx, iy);
            VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
            new_tristrip(&tlist, prev_edge, px, iy); 
            e1->outp[ABOVE]= prev_edge->outp[ABOVE];
            VERTEX(e1, ABOVE, RIGHT, ix, iy);
            new_tristrip(&tlist, e0, ix, iy);
            next_edge->outp[ABOVE]= e0->outp[ABOVE];
            VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
            break;
          case ILI:
            if (p)
            {
              VERTEX(e0, ABOVE, LEFT, ix, iy);
              N_EDGE(next_edge, e1, ABOVE, nx, iy);
              VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
              e1->outp[ABOVE]= e0->outp[ABOVE];
              e0->outp[ABOVE]= NULL;
            }
            break;
          case IRI:
            if (q)
            {
              VERTEX(e1, ABOVE, RIGHT, ix, iy);
              P_EDGE(prev_edge, e0, ABOVE, px, iy);
              VERTEX(prev_edge, ABOVE, LEFT, px, iy);
              e0->outp[ABOVE]= e1->outp[ABOVE];
              e1->outp[ABOVE]= NULL;
            }
            break;
          case IMX:
            if (p && q)
            {
              VERTEX(e0, ABOVE, RIGHT, ix, iy);
              VERTEX(e1, ABOVE, LEFT, ix, iy);
              e0->outp[ABOVE]= NULL;
              e1->outp[ABOVE]= NULL;
              P_EDGE(prev_edge, e0, ABOVE, px, iy);
              VERTEX(prev_edge, ABOVE, LEFT, px, iy);
              new_tristrip(&tlist, prev_edge, px, iy);
              N_EDGE(next_edge, e1, ABOVE, nx, iy);
              VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
              next_edge->outp[ABOVE]= prev_edge->outp[ABOVE];
              VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
            }
            break;
          case IMM:
            if (p && q)
            {
              VERTEX(e0, ABOVE, RIGHT, ix, iy);
              VERTEX(e1, ABOVE, LEFT, ix, iy);
              P_EDGE(prev_edge, e0, ABOVE, px, iy);
              VERTEX(prev_edge, ABOVE, LEFT, px, iy);
              new_tristrip(&tlist, prev_edge, px, iy);
              N_EDGE(next_edge, e1, ABOVE, nx, iy);
              VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
              e1->outp[ABOVE]= prev_edge->outp[ABOVE];
              VERTEX(e1, ABOVE, RIGHT, ix, iy);
              new_tristrip(&tlist, e0, ix, iy);
              next_edge->outp[ABOVE]= e0->outp[ABOVE];
              VERTEX(next_edge, ABOVE, RIGHT, nx, iy);
            }
            break;
          case EMM:
            if (p && q)
            {
              VERTEX(e0, ABOVE, LEFT, ix, iy);
              new_tristrip(&tlist, e1, ix, iy);
              e0->outp[ABOVE]= e1->outp[ABOVE];
            }
            break;
          default:
            break;
          } /* End of switch */
        } /* End of contributing intersection conditional */

        /* Swap bundle sides in response to edge crossing */
        if (e0->bundle[ABOVE][CLIP])
	        e1->bside[CLIP]= !e1->bside[CLIP];
        if (e1->bundle[ABOVE][CLIP])
	        e0->bside[CLIP]= !e0->bside[CLIP];
        if (e0->bundle[ABOVE][SUBJ])
	        e1->bside[SUBJ]= !e1->bside[SUBJ];
        if (e1->bundle[ABOVE][SUBJ])
	        e0->bside[SUBJ]= !e0->bside[SUBJ];

        /* Swap e0 and e1 bundles in the AET */
        prev_edge= e0->prev;
        next_edge= e1->next;
        if (e1->next)
          e1->next->prev= e0;

        if (e0->bstate[ABOVE] == BUNDLE_HEAD)
        {
          search= TRUE;
          while (search)
          {
            prev_edge= prev_edge->prev;
            if (prev_edge)
            {
              if (prev_edge->bundle[ABOVE][CLIP]
               || prev_edge->bundle[ABOVE][SUBJ]
               || (prev_edge->bstate[ABOVE] == BUNDLE_HEAD))
                search= FALSE;
            }
            else
              search= FALSE;
          }
        }
        if (!prev_edge)
        {
           e1->next= aet;
           aet= e0->next;
        }
        else
        {
          e1->next= prev_edge->next;
          prev_edge->next= e0->next;
        }
        e0->next->prev= prev_edge;
        e1->next->prev= e1;
        e0->next= next_edge;
      } /* End of IT loop*/

      /* Prepare for next scanbeam */
      for (edge= aet; edge; edge= next_edge)
      {
        next_edge= edge->next;
        succ_edge= edge->succ;

        if ((edge->top.y == yt) && succ_edge)
        {
          /* Replace AET edge by its successor */
          succ_edge->outp[BELOW]= edge->outp[ABOVE];
          succ_edge->bstate[BELOW]= edge->bstate[ABOVE];
          succ_edge->bundle[BELOW][CLIP]= edge->bundle[ABOVE][CLIP];
          succ_edge->bundle[BELOW][SUBJ]= edge->bundle[ABOVE][SUBJ];
          prev_edge= edge->prev;
          if (prev_edge)
            prev_edge->next= succ_edge;
          else
            aet= succ_edge;
          if (next_edge)
            next_edge->prev= succ_edge;
          succ_edge->prev= prev_edge;
          succ_edge->next= next_edge;
        }
        else
        {
          /* Update this edge */
          edge->outp[BELOW]= edge->outp[ABOVE];
          edge->bstate[BELOW]= edge->bstate[ABOVE];
          edge->bundle[BELOW][CLIP]= edge->bundle[ABOVE][CLIP];
          edge->bundle[BELOW][SUBJ]= edge->bundle[ABOVE][SUBJ];
          edge->xb= edge->xt;
        }
        edge->outp[ABOVE]= NULL;
      }
    }
  } /* === END OF SCANBEAM PROCESSING ================================== */

  /* Generate result tristrip from tlist */
  result->strip= NULL;
  result->num_strips= count_tristrips(tlist);
  if (result->num_strips > 0)
  {
    MY_NEW(result->strip, result->num_strips * sizeof(gpc_vertex_list),
           "tristrip list creation",gpc_vertex_list);

    s= 0;
    for (tn= tlist; tn; tn= tnn)
    {
      tnn= tn->next;

      if (tn->active > 2)
      {
        /* Valid tristrip: copy the vertices and free the heap */
        result->strip[s].num_vertices= tn->active;
        MY_NEW(result->strip[s].vertex, tn->active * sizeof(gpc_vertex),
               "tristrip creation",gpc_vertex);
        v= 0;
        if (INVERT_TRISTRIPS)
        {
          lt= tn->v[RIGHT];
          rt= tn->v[LEFT];
        }
        else
        {
          lt= tn->v[LEFT];
          rt= tn->v[RIGHT];
        }
        while (lt || rt)
        {
          if (lt)
          {
            ltn= lt->next;
            result->strip[s].vertex[v].x= lt->x;
            result->strip[s].vertex[v].y= lt->y;
            v++;
            MY_FREE(lt);
            lt= ltn;
          }
          if (rt)
          {
            rtn= rt->next;
            result->strip[s].vertex[v].x= rt->x;
            result->strip[s].vertex[v].y= rt->y;
            v++;
            MY_FREE(rt);
            rt= rtn;
          }
        }
        s++;
      }
      else
      {
        /* Invalid tristrip: just free the heap */
        for (lt= tn->v[LEFT]; lt; lt= ltn)
        {
          ltn= lt->next;
          MY_FREE(lt);
        }
        for (rt= tn->v[RIGHT]; rt; rt=rtn)
        {
          rtn= rt->next;
          MY_FREE(rt);
        }
      }
      MY_FREE(tn);
    }
  }

  /* Tidy up */
  reset_it(&it);
  reset_lmt(&lmt);
  MY_FREE(c_heap);
  MY_FREE(s_heap);
  MY_FREE(sbt);
}

/*
===========================================================================
                           End of file: gpc.c
===========================================================================
*/

GClipPolyBy_Vatti_algorithm::GClipPolyBy_Vatti_algorithm()
{
 ;
}

GClipPolyBy_Vatti_algorithm::~GClipPolyBy_Vatti_algorithm()
{
  ResetAll();
}

void GClipPolyBy_Vatti_algorithm::ResetAll             ()
{
  ResetSubjectPolys();
  ResetClipPolys   ();
  ResetResultPolys ();
}

void GClipPolyBy_Vatti_algorithm::ResetSubjectPolys    ()
{
  gpc_free_polygon(&m_SPolys);
}

void GClipPolyBy_Vatti_algorithm::ResetClipPolys       ()
{
  gpc_free_polygon(&m_CPolys);
}

void GClipPolyBy_Vatti_algorithm::ResetResultPolys     ()
{
  gpc_free_polygon(&m_RPolys);
}

//---------------------------------------------------------------------------------------------
// Clipping 대상 Polygon을 추가 한다. 
void GClipPolyBy_Vatti_algorithm::AddSubjectPoly (CArray<GPoint3D, GPoint3D&>& Poly,BOOL bIsHole)
{
  int IsHole = bIsHole ? 1 : 0;
  //int IsHole = bIsHole ? 0 : 1;
  gpc_vertex_list vlist;
  
  vlist.num_vertices = Poly.GetSize();

  MY_NEW(vlist.vertex,vlist.num_vertices*sizeof(gpc_vertex), "vertex creation",gpc_vertex);
  
  for(int i = 0 ; i < vlist.num_vertices; i++)
  {
    vlist.vertex[i].x = Poly[i].x;
    vlist.vertex[i].y = Poly[i].y;
  }

  gpc_add_contour(&m_SPolys,&vlist,IsHole);

  MY_FREE(vlist.vertex);
}

void GClipPolyBy_Vatti_algorithm::AddSubjectPoly (CArray<GPoint3D, GPoint3D>& Poly,BOOL bIsHole)
{
  int IsHole = bIsHole ? 1 : 0;

  //int IsHole = bIsHole ? 0 : 1;
  gpc_vertex_list vlist;
  
  vlist.num_vertices = Poly.GetSize();

  MY_NEW(vlist.vertex,vlist.num_vertices*sizeof(gpc_vertex), "vertex creation",gpc_vertex);
  
  for(int i = 0 ; i < vlist.num_vertices; i++)
  {
    vlist.vertex[i].x = Poly[i].x;
    vlist.vertex[i].y = Poly[i].y;
  }

  gpc_add_contour(&m_SPolys,&vlist,IsHole);

  MY_FREE(vlist.vertex);
}

void GClipPolyBy_Vatti_algorithm::AddSubjectPoly       (CArray<DGL_3dp,DGL_3dp&>& Poly,BOOL bIsHole)
{
  int IsHole = bIsHole ? 1 : 0;

  //int IsHole = bIsHole ? 0 : 1;
  gpc_vertex_list vlist;
  
  vlist.num_vertices = Poly.GetSize();

  MY_NEW(vlist.vertex,vlist.num_vertices*sizeof(gpc_vertex), "vertex creation",gpc_vertex);
  
  for(int i = 0 ; i < vlist.num_vertices; i++)
  {
    vlist.vertex[i].x = Poly[i].x();
    vlist.vertex[i].y = Poly[i].y();
  }

  gpc_add_contour(&m_SPolys,&vlist,IsHole);

  MY_FREE(vlist.vertex);
}


//----------------------------------------------------------------------------------------------
// Clipping 영역을 지정한다.
void GClipPolyBy_Vatti_algorithm::AddClipPoly    (CArray<GPoint3D, GPoint3D&>& Poly)
{
  gpc_vertex_list vlist;
  
  vlist.num_vertices = Poly.GetSize();

  MY_NEW(vlist.vertex,vlist.num_vertices*sizeof(gpc_vertex), "vertex creation",gpc_vertex);
  
  for(int i = 0 ; i < vlist.num_vertices; i++)
  {
    vlist.vertex[i].x = Poly[i].x;
    vlist.vertex[i].y = Poly[i].y;
  }

  gpc_add_contour(&m_CPolys,&vlist,0);

  MY_FREE(vlist.vertex);
}

void GClipPolyBy_Vatti_algorithm::AddClipPoly    (CArray<GPoint3D, GPoint3D>& Poly)
{
  gpc_vertex_list vlist;
  
  vlist.num_vertices = Poly.GetSize();

  MY_NEW(vlist.vertex,vlist.num_vertices*sizeof(gpc_vertex), "vertex creation",gpc_vertex);
  
  for(int i = 0 ; i < vlist.num_vertices; i++)
  {
    vlist.vertex[i].x = Poly[i].x;
    vlist.vertex[i].y = Poly[i].y;
  }

  gpc_add_contour(&m_CPolys,&vlist,0);

  MY_FREE(vlist.vertex);
}

void  GClipPolyBy_Vatti_algorithm::AddClipPoly          (CArray<DGL_3dp,DGL_3dp&>& Poly)
{
  gpc_vertex_list vlist;
  
  vlist.num_vertices = Poly.GetSize();

  MY_NEW(vlist.vertex,vlist.num_vertices*sizeof(gpc_vertex), "vertex creation",gpc_vertex);
  
  for(int i = 0 ; i < vlist.num_vertices; i++)
  {
    vlist.vertex[i].x = Poly[i].x();
    vlist.vertex[i].y = Poly[i].y();
  }

  gpc_add_contour(&m_CPolys,&vlist,0);

  MY_FREE(vlist.vertex);
}

void  GClipPolyBy_Vatti_algorithm::MakeResult_intersect()
{
  ResetResultPolys();
  gpc_polygon_clip(GPC_INT,
                   &m_SPolys,
                   &m_CPolys,
                   &m_RPolys);
}

void  GClipPolyBy_Vatti_algorithm::MakeResult_merge()
{
  ResetResultPolys();
  gpc_polygon_clip(GPC_UNION,
                   &m_SPolys,
                   &m_CPolys,
                   &m_RPolys);
}

void GClipPolyBy_Vatti_algorithm::MakeResult_extract   ()
{
  ResetResultPolys();
  gpc_polygon_clip(GPC_DIFF ,
                   &m_SPolys,
                   &m_CPolys,
                   &m_RPolys);
}

int  GClipPolyBy_Vatti_algorithm::GetResultPolyIndex(BOOL bIsHole,CArray<int,int>&arIPolys)
{
  arIPolys.RemoveAll();
  int nP = m_RPolys.num_contours;
  for(int i = 0; i < nP; i++)
  {
    if(bIsHole == TRUE && m_RPolys.hole[i] == 1)
    {
      arIPolys.Add(i);
    }
    else if(bIsHole == FALSE && m_RPolys.hole[i] == 0)
    {
      arIPolys.Add(i);
    }
    else
    {
      //int i = 0; 
      //i++;
    }
  }
  return arIPolys.GetSize();;
}

BOOL GClipPolyBy_Vatti_algorithm::GetResultPoly(int nIndex,CArray<GPoint3D, GPoint3D&>&Poly)
{
  if(m_RPolys.num_contours -1 < nIndex || nIndex < 0) return FALSE;
  Poly.RemoveAll();
 
  gpc_vertex_list* pVL;
  pVL = m_RPolys.contour+nIndex;
  int nV = pVL->num_vertices;
  
  GPoint3D Pos;
  Poly.SetSize(nV);
  for(int i = 0 ; i < nV; i++)
  {
    Pos.Set(pVL->vertex[i].x,pVL->vertex[i].y,0);
    Poly.SetAt(i,Pos);
  }
  return TRUE;
}

BOOL GClipPolyBy_Vatti_algorithm::GetResultPoly(int nIndex,CArray<GPoint3D, GPoint3D>&Poly)
{
  if(m_RPolys.num_contours -1 < nIndex || nIndex < 0) return FALSE;
  Poly.RemoveAll();
 
  gpc_vertex_list* pVL;
  pVL = m_RPolys.contour+nIndex;
  int nV = pVL->num_vertices;
  
  GPoint3D Pos;
  Poly.SetSize(nV);
  for(int i = 0 ; i < nV; i++)
  {
    Pos.Set(pVL->vertex[i].x,pVL->vertex[i].y,0);
    Poly.SetAt(i,Pos);
  }
  return TRUE;
}

BOOL GClipPolyBy_Vatti_algorithm::GetResultPoly     (int nIndex,CArray<DGL_3dp,DGL_3dp&>& Poly)
{
  if(m_RPolys.num_contours -1 < nIndex || nIndex < 0) return FALSE;
  Poly.RemoveAll();
 
  gpc_vertex_list* pVL;
  pVL = m_RPolys.contour+nIndex;
  int nV = pVL->num_vertices;
  
  DGL_3dp Pos;
  Poly.SetSize(nV);
  for(int i = 0 ; i < nV; i++)
  {
    Pos.Set(pVL->vertex[i].x,pVL->vertex[i].y,0);
    Poly.SetAt(i,Pos);
  }
  return TRUE;
}



// GPoint3D * arg1, GPoint3D * arg2
int Vatti__CompareXofPoint( const void* arg1, const void * arg2)
{
  GPoint3D P1,P2;
  P1 = *(GPoint3D*)arg1;
  P2 = *(GPoint3D*)arg2;

  if(GMU::IsSameDoubleV(P1.x,P2.x)) return 0;
  if(P1.x < P2.x) return -1;
  return 1;
}

int Vatti__CompareYofPoint( const void* arg1, const void * arg2)
{
  GPoint3D P1,P2;
  P1 = *(GPoint3D*)arg1;
  P2 = *(GPoint3D*)arg2;

  if(GMU::IsSameDoubleV(P1.y,P2.y)) return 0;
  if(P1.y < P2.y) return -1;
  return 1;
}

int Vatti__CompareZofPoint( const void* arg1, const void * arg2)
{
  GPoint3D P1,P2;
  P1 = *(GPoint3D*)arg1;
  P2 = *(GPoint3D*)arg2;

  if(GMU::IsSameDoubleV(P1.z,P2.z)) return 0;
  if(P1.z < P2.z) return -1;
  return 1;
}

//----------------------------------------------------------------------------
// Vatti Algolithm 구현에서는 Cutting Edge를 구할 수 없다. 
// 별도 구현한다. 

bool GClipPolyBy_Vatti_algorithm::DistFromLineToPoint(GPoint3D &LP1,GPoint3D &LP2,GPoint3D &FromP,double & Dist)
{
  
  if(LP1 == LP2) 
  {
    double Dx   = LP1.x - FromP.x;
    double Dy   = LP1.y - FromP.y;
    Dist = sqrt(Dx*Dx + Dy*Dy);
    return false;
  }
  
  double 	a,b,c,a1,b1,c1;
	double	x,y,x1,y1,x2,y2;
	double 	Test;//,Test2;
	double	Temp_3;
  GPoint3D P1 = LP1;
  GPoint3D P2 = LP2;
	
	x = FromP.x;
	y = FromP.y;

	x1 = P1.x;
	y1 = P1.y;
	
	x2 = P2.x;
	y2 = P2.y;

	a=y2-y1;
	b=x1-x2;
	c=x2*y1-x1*y2;
	
	a1=b;
	b1=-a;
	c1=a*y-b*x;
	
	Test  =(a1*x1+b1*y1+c1)*(a1*x2+b1*y2+c1);
  //Test2 =(a*x+b*y+c);
  //if(GMU::IsZero(Test2))
  //{
  //  Dist = 0;   
  //  return TRUE;
  //}
	
	if (Test<0)
	{	
		Temp_3=(double)(a*x+b*y+c);
    Temp_3=Temp_3*Temp_3;
		Dist =sqrt(Temp_3/(a*a+b*b));
    return true;
	}
  
  return false;
}


bool GClipPolyBy_Vatti_algorithm::GetLineIntersectPoint2D(GLine3D &LToTest, GLine3D& TempL, GPoint3D &IntersectPoint)
{
	GPoint3D P1,P2,P3,P4;

  P1 = LToTest.P1; P2 = LToTest.P2; 
  P3 = TempL.P1  ; P4 = TempL.P2  ;
    
  double A,B,C,D;
	double T,V;
	double DivT,DivV;
	
	A = P4.x - P3.x;
	B = P2.x - P1.x;
	C = P4.y - P3.y;
	D = P2.y - P1.y;

	DivT = A*D - B*C; 
	DivV = B*C - A*D;
	
  if(!GMU::IsZero(DivT) && !GMU::IsZero(DivV))
	{
		T = (A*P3.y + C*P1.x - C*P3.x - A*P1.y) / DivT;
		V = (B*P1.y + D*P3.x - B*P3.y - D*P1.x) / DivV;

		//if( (T > -V_LGL_LINE_INTERSECT_EPSILON && T < 1.0 + V_LGL_LINE_INTERSECT_EPSILON) &&
		//	(V > -V_LGL_LINE_INTERSECT_EPSILON && V < 1.0 + V_LGL_LINE_INTERSECT_EPSILON))
		{
			IntersectPoint.x = P1.x + B * T;
			IntersectPoint.y = P1.y + D * T;
			IntersectPoint.z = 0.0;
	
			return true;
		}
	}

  //-------------------------------------------------------------------------------------------
	IntersectPoint.Set(0,0,0);
	return false;
}

//---------------------------------------------------------------------------------------------
// Edge에 걸린것은 포함되지 않은것으로 처리한다. 
bool GClipPolyBy_Vatti_algorithm::IsWithin2DByAngle(gpc_vertex_list & aPolygon, GPoint3D point,bool& bIsOnEdge)
{
  
  GPoint3D MinP,MaxP;
  if(aPolygon.num_vertices < 3) return false;
  MinP.Set(aPolygon.vertex[0].x,aPolygon.vertex[0].y,0);
  MaxP.Set(aPolygon.vertex[0].x,aPolygon.vertex[0].y,0);
  
  GPoint3D TempP;
  int VertexNum = aPolygon.num_vertices;

  for(int i = 1; i < VertexNum; i++)
  {
		TempP.Set(aPolygon.vertex[i].x,aPolygon.vertex[i].y,0);
		
		if(MinP.x > TempP.x) MinP.x = TempP.x;
		if(MinP.y > TempP.y) MinP.y = TempP.y;
		if(MaxP.x < TempP.x) MaxP.x = TempP.x;
		if(MaxP.y < TempP.y) MaxP.y = TempP.y;
  }

  GGeometryEngine GM;

  GRect OR;
  OR.Set(MinP.x,MinP.y,MaxP.x,MaxP.y);
  //----------------------------------------------------------------
  //폴리곤에 완전히 포함 되지 않는 point는 포함 테스트를 생략한다.
  if(! GM._IsWithin2D(OR,&point)) return FALSE;

  CArray<GVector, GVector&> Vectors;
  GVector TVector;

  //----------------------------------------------------------------
  // Vertext또는 Edge에 포함된 점 관련 처리 //
  
  for(int i = 0 ; i < VertexNum ; i++)
  {
    TempP.Set(aPolygon.vertex[i].x,aPolygon.vertex[i].y,0);
    TVector.Set(point,TempP);
    TVector.MakeUnit();
    
    Vectors.Add(TVector);
  }

  int nVectors = Vectors.GetSize();

  double AngleSum, TAngle;
  AngleSum = 0;
  GVector NVector;

  double ValPI  = M_PI;

  bIsOnEdge = false;

  for(int i = 0 ; i < nVectors ; i++)
  {
    
    if(i == nVectors-1)
    {
      if(Vectors[i].Angle3D(Vectors[0],TAngle))
      {
        if(!GMU::IsZero(TAngle))
        {
          if(!GMU::IsSameDouble(TAngle,ValPI)) // Edge에 걸린경우는 포함되지 않은것으로 처리한다. 
          {
            NVector = Vectors[i].Cross(Vectors[0]);
            if(NVector.z < 0) TAngle = -TAngle;
            AngleSum += TAngle;
          }
          else
          {
            bIsOnEdge = true;
          }
        }
      }   
    }
    else
    {
      if(Vectors[i].Angle3D(Vectors[i+1],TAngle))
      {
        if(!GMU::IsZero(TAngle))
        {
          if(!GMU::IsSameDouble(TAngle,ValPI)) // Edge에 걸린경우는 포함되지 않은것으로 처리한다. 
          {
            NVector = Vectors[i].Cross(Vectors[i+1]);
            if(NVector.z < 0) TAngle = -TAngle;
            AngleSum += TAngle;
          }
          else
          {
            bIsOnEdge = true;
          }
        }
      }
    }
  }
  
  double Val2PI  = M_PI *  2;
  double MVal2PI = M_PI * -2;

  AngleSum = fabs(AngleSum);
  if(fabs(AngleSum - Val2PI) <= 0.0001) return true;

  return false;
}

bool GClipPolyBy_Vatti_algorithm::GetCuttingEdge(CArray<CuttingEdge,CuttingEdge&> &Edges,GLine3D& CuttingLine)
{
  CArray<gpc_vertex_list* , gpc_vertex_list*> arOPolys;
  CArray<gpc_vertex_list* , gpc_vertex_list*> arIPolys;
  CArray<gpc_vertex_list* , gpc_vertex_list*> arAllPolys;
  int nPolys = m_SPolys.num_contours; // Outter , Holes ...
  
  if(nPolys == 0) return false;

  Edges.RemoveAll();

  for(int i = 0; i < nPolys ; i++)
  {
    // Polygon 조건을 만족하지 않으면 오류 처리한다. 
    if(m_SPolys.contour[i].num_vertices < 3) 
      return false; 
    
    if(m_SPolys.hole[i])
      arIPolys.Add(&(m_SPolys.contour[i]));
    else
      arOPolys.Add(&(m_SPolys.contour[i]));

    arAllPolys.Add(&(m_SPolys.contour[i]));
  }

  gpc_vertex_list *pCurPoly;
  GPoint3D MinP,MaxP,TempP,IntersectPoint;
  MinP.x = MinP.y = DBL_MAX ;
  MaxP.x = MaxP.y = -DBL_MAX;
  MinP.z = 0.0;
  MaxP.z = 0.0;
  int VertexNum;
  
  for(int i = 0; i < arIPolys.GetSize(); i ++)
  {
    pCurPoly = arIPolys[i];
    VertexNum = pCurPoly->num_vertices;
    
    for(int j = 0; j < VertexNum; j++)
    {
			TempP.Set(pCurPoly->vertex[j].x, pCurPoly->vertex[j].y,0);
      
      MinP.x = min(MinP.x,TempP.x);
      MinP.y = min(MinP.y,TempP.y);
      MaxP.x = max(MaxP.x,TempP.x);
      MaxP.y = max(MaxP.y,TempP.y);
    }
  }

  GRect OR;
  OR.Set(MinP.x,MinP.y,MaxP.x,MaxP.y);
  //Cutting Line의 끝점이 항상 최외곽 폴리곤의 외부에 있도록 임의 조정할것. 

  CArray<GPoint3D,GPoint3D> IntsPoints;
  
  IntsPoints.Add(CuttingLine.P1);

  int nAllPolys = arAllPolys.GetSize();
  int IP1,IP2; 
  GLine3D TempL;
  
  //설정된 모든 Polygon과 Cutting Line과의 교점을 구한다. 
  for(int i = 0; i < nAllPolys ; i++)
  {
    pCurPoly = arAllPolys[i];
    VertexNum = pCurPoly->num_vertices;

    for(int j = 0; j < VertexNum; j++)
    {
      if(j == VertexNum -1)
      {
        IP1 = j ; IP2 = 0;
      } 
  	  else
      {
        IP1 = j ; IP2 = j+1;
      }
      
      TempL.P1.Set(pCurPoly->vertex[IP1].x,pCurPoly->vertex[IP1].y,0);
      TempL.P2.Set(pCurPoly->vertex[IP2].x,pCurPoly->vertex[IP2].y,0);
      if(GetLineIntersectPoint2D(CuttingLine,TempL,IntersectPoint))
      {
        IntsPoints.Add(IntersectPoint);
      }
    }
  }

  IntsPoints.Add(CuttingLine.P2);

  //---------------------------------------------------------------------------------------------------
  // 교차점들로 생성된 Line들의 중점이 대상 Polygon에 포함되는지 검사한다. 
  // 중점이 Polygon 내부에 존재하면 그 Line은 명벽히 Polygon내부에 포함된Line이다 
  // 즉 Clipping 처리된 결과이다.
  {
    int nInts  = IntsPoints.GetSize();
    int nIPoly = arIPolys.GetSize();
    int nOPoly = arOPolys.GetSize();
    
    if(nInts >= 2)
    {
      // XSorting
      double Dx,Dy;
      GPoint3D ICentPoint;
      gpc_vertex_list *pPolygon;
      CuttingEdge CutEdge;

      Dx =fabs(CuttingLine.P1.x - CuttingLine.P2.x);
      Dy =fabs(CuttingLine.P1.y - CuttingLine.P2.y);
      
      if(Dy > Dx)
        qsort(IntsPoints.GetData(),nInts,sizeof(GPoint3D),Vatti__CompareYofPoint);
      else
        qsort(IntsPoints.GetData(),nInts,sizeof(GPoint3D),Vatti__CompareXofPoint);
  
      bool bIsOnEdge_I, bIsOnEdge_O;
      for(int i = 0 ; i < nInts-1; i++)
      {
        if(IntsPoints[i] == IntsPoints[i+1]) continue;
        
        ICentPoint = IntsPoints[i] + IntsPoints[i+1];
                  
        ICentPoint.x = ICentPoint.x /2.;
        ICentPoint.y = ICentPoint.y /2.;
        ICentPoint.z = 0;
       
        
        for(int j = 0; j < nOPoly; j++)
        {
          pPolygon = arOPolys[j];
          bIsOnEdge_O = bIsOnEdge_I = false;
          if(IsWithin2DByAngle(*pPolygon,ICentPoint,bIsOnEdge_O))
          {
            if(bIsOnEdge_O) break;
            
            bool bIsInnerIPoly = false;
            for(int k = 0; k < nIPoly; k++)
            {
              pPolygon = arIPolys[k];
              if(!IsWithin2DByAngle(*pPolygon, ICentPoint,bIsOnEdge_I))
              {
                if(bIsOnEdge_I) //Hole의 Edge에 존재하는 Point는 제외시킨다. 
                {
                  bIsInnerIPoly = true;
                }
              }
              else
              {
                bIsInnerIPoly = true;
              }
            }

            if(false == bIsInnerIPoly)
            {
              CutEdge.m_Pos1 = IntsPoints[i];
              CutEdge.m_Pos2 = IntsPoints[i+1];
              if(CutEdge.IsValid())
                Edges.Add(CutEdge);
            }
          }
        }
      }
    }
  }

  return true;
}


/*

//-----------------------------------------------------------------------------
//  
// 
VattisPolyClipEventHandler::VattisPolyClipEventHandler()
{
  m_nMode = -1;
  m_pCurPolyV = NULL;
  SetOperationIntersect();
}

VattisPolyClipEventHandler::~VattisPolyClipEventHandler()
{
  ResetAllData();
}

void VattisPolyClipEventHandler::StartSPolyInput(BOOL bIsHole)
{
  m_nMode = -1;
  
  if(!bIsHole)
  {
    int nPoly = m_arOPolys.GetSize();
    if(nPoly != 0)
    {
      int nVert = m_arOPolys[nPoly-1]->GetSize();
      if(nVert < 3)
        m_arOPolys[nPoly-1]->RemoveAll();
      else
      {
        if(m_arOPolys[nPoly-1]->GetAt(nVert-1) != 
           m_arOPolys[nPoly-1]->GetAt(0))
        {
          m_arOPolys[nPoly-1]->RemoveAll();
        }
        else
        {
          m_pCurPolyV = new CArray<GPoint3D,GPoint3D>;
          m_arOPolys.Add(m_pCurPolyV);
        }
      }
    }
    else
    {
      m_pCurPolyV = new CArray<GPoint3D,GPoint3D>;
      m_arOPolys.Add(m_pCurPolyV);
    }
    m_nMode = 0;
  }
  else
  {
    int nPoly = m_arIPolys.GetSize();
    if(nPoly != 0)
    {
      int nVert = m_arIPolys[nPoly-1]->GetSize();
      if(nVert < 3)
        m_arIPolys[nPoly-1]->RemoveAll();
      else
      {
        if(m_arIPolys[nPoly-1]->GetAt(nVert-1) != 
           m_arIPolys[nPoly-1]->GetAt(0))
        {
          m_arIPolys[nPoly-1]->RemoveAll();
        }
        else
        {
          m_pCurPolyV = new CArray<GPoint3D,GPoint3D>;
          m_arIPolys.Add(m_pCurPolyV);
        }
      }
    }
    else
    {
      m_pCurPolyV = new CArray<GPoint3D,GPoint3D>;
      m_arIPolys.Add(m_pCurPolyV);  
    }
    m_nMode = 1;    
  }
}

void VattisPolyClipEventHandler::StartCPolyInput()
{
  m_nMode = -1;
  int nPoly = m_arCPolys.GetSize();
  if(nPoly != 0)
  {
    int nVert = m_arCPolys[nPoly-1]->GetSize();
    if(nVert < 3)
      m_arCPolys[nPoly-1]->RemoveAll();
    else
    {
      if(m_arCPolys[nPoly-1]->GetAt(nVert-1) != 
         m_arCPolys[nPoly-1]->GetAt(0))
      {
        m_arCPolys[nPoly-1]->RemoveAll();
      }
      else
      {
        m_pCurPolyV = new CArray<GPoint3D,GPoint3D>;
        m_arCPolys.Add(m_pCurPolyV);
      }
    }
  }
  else
  {
    m_pCurPolyV = new CArray<GPoint3D,GPoint3D>;
    m_arCPolys.Add(m_pCurPolyV);
  }
  m_nMode = 2;
}

void VattisPolyClipEventHandler::RemoveCPolys()
{
  m_nMode = -1;
  int nPolys = m_arCPolys.GetSize();
  for(int i = 0; i < nPolys; i++)
  {
    delete m_arCPolys[i];
  }
  m_arCPolys.RemoveAll();
  m_pCurPolyV = NULL;
}

void VattisPolyClipEventHandler::RemoveSPolys(BOOL bIsHole)
{
  m_nMode = -1;
  int nPolys;
  if(bIsHole == FALSE)
  {
    nPolys = m_arOPolys.GetSize();
    for(int i = 0; i < nPolys; i++)
    {
      delete m_arOPolys[i];
    }
    m_arOPolys.RemoveAll();
  }

  if(bIsHole == TRUE)
  {
    nPolys = m_arIPolys.GetSize();
    for(int i = 0; i < nPolys; i++)
    {
      delete m_arIPolys[i];
    }
    m_arIPolys.RemoveAll();
  }
  m_pCurPolyV = NULL;
}

void VattisPolyClipEventHandler::RemoveRPolys()
{
  m_nMode = -1;
  int nPolys = m_arRPolys.GetSize();
  for(int i = 0; i < nPolys; i++)
  {
    delete m_arRPolys[i];
  }
  m_arRPolys.RemoveAll();

  nPolys = m_arHPolys.GetSize();
  for(    i = 0; i < nPolys; i++)
  {
    delete m_arHPolys[i];
  }
  m_arHPolys.RemoveAll();
  m_pCurPolyV = NULL;
}

void VattisPolyClipEventHandler::ResetAllData()
{
  RemoveSPolys(TRUE ); // Hole 
  RemoveSPolys(FALSE); // Normal
  RemoveCPolys();
  RemoveRPolys();
}

//----------------------------------------------------------------  
void VattisPolyClipEventHandler::MakeResult(GViewWnd* pVW)
{
  RemoveRPolys();
  
  GClipPolyBy_Vatti_algorithm Vatti;
  
  int nPoly = m_arIPolys.GetSize();
  for(int i = 0 ; i < nPoly; i++)
  {
    Vatti.AddSubjectPoly(*m_arIPolys[i],TRUE);
  }

  nPoly = m_arOPolys.GetSize();
  for(    i = 0; i < nPoly; i++)
  {
    Vatti.AddSubjectPoly(*m_arOPolys[i],FALSE);
  }

  nPoly = m_arCPolys.GetSize();
  for(    i = 0; i < nPoly; i++)
  {
    Vatti.AddClipPoly(*m_arCPolys[i]);
  }
  
  if(IsOperationIntersect())
    Vatti.MakeResult_intersect();
  else if(IsOperationMerge())
    Vatti.MakeResult_merge();
  else if(IsOperationExtract())
    Vatti.MakeResult_extract();
  else 
  {
    ASSERT(FALSE);
  }

  CArray<int,int>arPolyI;
  CArray<GPoint3D, GPoint3D>* pPoly;
  Vatti.GetResultPolyIndex(TRUE ,arPolyI);
  nPoly = arPolyI.GetSize();
  for(    i = 0; i < nPoly; i++)
  {
    pPoly = new CArray<GPoint3D, GPoint3D>;
    Vatti.GetResultPoly(arPolyI[i],*pPoly);
    m_arHPolys.Add(pPoly);
    if(pPoly->GetSize())
      pPoly->Add((*pPoly)[0]);
  }

  arPolyI.RemoveAll();
  Vatti.GetResultPolyIndex(FALSE,arPolyI);
  nPoly = arPolyI.GetSize();
  for(    i = 0; i < nPoly; i++)
  {
    pPoly = new CArray<GPoint3D, GPoint3D>;
    Vatti.GetResultPoly(arPolyI[i],*pPoly);
    m_arRPolys.Add(pPoly);
    if(pPoly->GetSize())
      pPoly->Add((*pPoly)[0]);
  }
  
  if(m_arCPolys.GetSize() && (m_arCPolys[0]->GetSize()>=3))
  {
    GLine3D CuttingLine;
    CuttingLine.Set((*m_arCPolys[0])[0],(*m_arCPolys[0])[1]);
    Vatti.GetCuttingEdge(m_arCuttingEdge,CuttingLine);
  }
  else
  {
    m_arCuttingEdge.RemoveAll();
  }
  
  pVW->Invalidate(); 
}


BOOL VattisPolyClipEventHandler::IsSamePoint(GPoint3D* pP1,GPoint3D* pP2, int Tol,GViewWnd* pVW)
{
  GLine3D Line;
  GLine2D CLine;
  Line.P1 = *pP1; Line.P2 = *pP2; 
  pVW->GetRE()->GetLineOnClientWindow(Line,CLine);
  if(Tol >= CLine.Length()) return TRUE;
  return FALSE;
}


bool VattisPolyClipEventHandler::GetBoundingBox(GPoint3D& MinP, GPoint3D& MaxP)
{

  int nOPolys = m_arOPolys.GetSize();
  int nIPolys = m_arIPolys.GetSize();

  if(nOPolys == 0) return false;

  double MinX, MinY, MaxX ,MaxY;

  MinX = MinY = DBL_MAX;
  MaxX = MaxY = -DBL_MAX;
 

  for(int i = 0 ; i < nOPolys; i++)
  {
    int nOVert = m_arOPolys[i]->GetSize();
    for( int j = 0; j < nOVert; j++)
    {
      MinX = min((*m_arOPolys[i])[j].x,MinX); 
      MinY = min((*m_arOPolys[i])[j].y,MinY); 
      MaxX = max((*m_arOPolys[i])[j].x,MaxX); 
      MaxY = max((*m_arOPolys[i])[j].y,MaxY); 
    }
  }

  for( i = 0; i < nIPolys; i++)
  {
    int nIVert = m_arIPolys[i]->GetSize();
    for( int j = 0; j < nIVert; j++)
    {
      MinX = min((*m_arIPolys[i])[j].x,MinX); 
      MinY = min((*m_arIPolys[i])[j].y,MinY); 
      MaxX = max((*m_arIPolys[i])[j].x,MaxX); 
      MaxY = max((*m_arIPolys[i])[j].y,MaxY);  
    }
  }

  MinP.Set(MinX,MinY,0);
  MaxP.Set(MaxX,MaxY,0);
  
  return TRUE;

}  
//----------------------------------------------------------------
//  Event Handler !!!
void VattisPolyClipEventHandler::On_LButtonDown(GPoint3D *pWcsPos,GViewWnd* pVW)
{
  int Tol = 2;
  
  //---------------------------------------------------------------------------------
  //(0) Input Outter Polygon (1) Input Inner Polygon (3) Input Clip Polygon  (4) Display Result 
  switch(m_nMode)
  {
  case 0: // Input Outter Polygon 
    {
      int nIPA = m_arOPolys.GetSize();
      if(nIPA)
      {
        CArray<GPoint3D,GPoint3D> *pAR;
        pAR = m_arOPolys[nIPA-1];
        int nV = pAR->GetSize();
        if(nV >= 3)
        {
          if((*pAR)[0] == (*pAR)[nV-1]) // 종료된것.
          {
            break;  
          }
        }

        if(nV >= 2)
        {
          if(IsSamePoint(&(*pAR)[0],pWcsPos,Tol,pVW))
            pAR->Add((*pAR)[0]);
          else
            pAR->Add(*pWcsPos);
        }
        else
        {
          pAR->Add(*pWcsPos);
        }

        pVW->Invalidate();
      }
    }
    break;
  case 1: // Input Inner Polygon 
    {
      int nIPA = m_arIPolys.GetSize();
      if(nIPA)
      {
        CArray<GPoint3D,GPoint3D> *pAR;
        pAR = m_arIPolys[nIPA-1];
        int nV = pAR->GetSize();
        if(nV >= 3)
        {
          if((*pAR)[0] == (*pAR)[nV-1]) // 종료된것.
          {
            break;  
          }
        }
        
        if(nV >= 2)
        {
          if(IsSamePoint(&(*pAR)[0],pWcsPos,Tol,pVW))
            pAR->Add((*pAR)[0]);
          else
            pAR->Add(*pWcsPos);
        }
        else
        {
          pAR->Add(*pWcsPos);
        }
        pVW->Invalidate();
      }
    }
    break;
  case 3: // Input Clip Polygon 
    {
      int nIPA = m_arCPolys.GetSize();
      if(nIPA)
      {
        CArray<GPoint3D,GPoint3D> *pAR;
        pAR = m_arCPolys[nIPA-1];
        int nV = pAR->GetSize();
        if(nV >= 3)
        {
          if((*pAR)[0] == (*pAR)[nV-1]) // 종료된것.
          {
            break;  
          }
        }
        
        if(nV >= 2)
        {
          if(IsSamePoint(&(*pAR)[0],pWcsPos,Tol,pVW))
            pAR->Add((*pAR)[0]);
          else
            pAR->Add(*pWcsPos);
        }
        else
        {
          pAR->Add(*pWcsPos);
        }
        pVW->Invalidate();
      }
    }
    break;
  case 2: // Input Clip Line
    {
      int nIPA = m_arCPolys.GetSize();
      if(nIPA)
      {
        CArray<GPoint3D,GPoint3D> *pAR;
        pAR = m_arCPolys[nIPA-1];
        int nV = pAR->GetSize();
        if(nV >= 3)
        {
          if((*pAR)[0] == (*pAR)[nV-1]) // 종료된것.
          {
            break;  
          }
        }
        
        if(nV >= 2)
        {
          
        }
        else
        {
          pAR->Add(*pWcsPos);
          if(pAR->GetSize() == 2)
          {
            GVector CLineVector, YVector,ZVector;
            GPoint3D P1, P2,P3,P4;

            P1.Set((*pAR)[0].x,(*pAR)[0].y,0);
            P2.Set(pWcsPos->x,pWcsPos->y,0);

            GPoint3D MinP,MaxP;
            GLine3D TempLine;
            GetBoundingBox(MinP,MaxP);

            TempLine.P1.Set(MinP.x,MinP.y,0);
            TempLine.P2.Set(MaxP.x,MaxP.y,0);

            double LOffset, YOffset;
            LOffset = TempLine.Length();
            YOffset = LOffset*2.;
            
            ZVector.Set(0,0,1);
            CLineVector.Set(P1,P2); CLineVector.MakeUnit();
            YVector = ZVector * CLineVector;

            CLineVector.ScalingMySelf(LOffset);

            P1.x = P1.x - CLineVector.x;
            P1.y = P1.y - CLineVector.y;
            P2.x = P2.x + CLineVector.x;
            P2.y = P2.y + CLineVector.y;
              
            YVector.ScalingMySelf(YOffset);

            P3.x = P2.x + YVector.x;
            P3.y = P2.y + YVector.y;
            P4.x = P1.x + YVector.x;
            P4.y = P1.y + YVector.y;
            
            pAR->RemoveAll();
            pAR->Add(P1);pAR->Add(P2);pAR->Add(P3);pAR->Add(P4);pAR->Add(P1);
          }
        }
        pVW->Invalidate();
      }
    }
    break;
  default:
    break;
  }
}

//----------------------------------------------------------------
// Paint Handler !!!
//
void VattisPolyClipEventHandler::On_Paint(GViewWnd* pVW,CDC*pDC)
{
  pDC->SaveDC();
  
  int nPoly,i;
  
  nPoly = this->m_arOPolys.GetSize();
  CPen OutPen(PS_SOLID,1,BLUE);
  pDC->SelectObject(&OutPen);
  for(    i = 0; i < nPoly; i++)
  {
    pVW->GetRE()->DrawPLine3DOnClientWindow(pDC,*m_arOPolys[i]);
  }
  
  nPoly = this->m_arIPolys.GetSize();
  CPen InnPen(PS_SOLID,1,RED);
  pDC->SelectObject(&InnPen);
  for(    i = 0; i < nPoly; i++)
  {
    pVW->GetRE()->DrawPLine3DOnClientWindow(pDC,*m_arIPolys[i]);
  }

  nPoly = this->m_arCPolys.GetSize();
  CPen CLPen(PS_SOLID,1,MEDIUM_GRAY);
  pDC->SelectObject(&CLPen);
  for(    i = 0; i < nPoly; i++)
  {
    pVW->GetRE()->DrawPLine3DOnClientWindow(pDC,*m_arCPolys[i]);
  }

  nPoly = this->m_arRPolys.GetSize();
  CPen   ROPen(PS_SOLID,2,BLACK);
  //CBrush ROBrush(WHITE);
  pDC->SelectObject(&ROPen);
  //pDC->SelectObject(&ROBrush);
  for(    i = 0; i < nPoly; i++)
  {
    pVW->GetRE()->DrawPLine3DOnClientWindow(pDC,*m_arRPolys[i]);
    //pVW->GetRE()->DrawPolygon3DFillOnClientWindow(pDC,*m_arRPolys[i]);
  }

  nPoly = this->m_arHPolys.GetSize();
  //CPen   RHPen  (PS_NULL,1,RED);
  CPen   RHPen  (PS_SOLID,1,RED);
  //CBrush RHBrush(WHITE);
  pDC->SelectObject(&RHPen);
  //pDC->SelectObject(&RHBrush);
  for(   i = 0; i < nPoly; i++)
  {
    pVW->GetRE()->DrawPLine3DOnClientWindow(pDC,*m_arHPolys[i]);
    //pVW->GetRE()->DrawPolygon3DFillOnClientWindow(pDC,*m_arHPolys[i]);
  }


  //-------------------------------------------------------------------
  // cutting Edge 표시 
  int nLine = m_arCuttingEdge.GetSize();
  CPen   CutLPen  (PS_SOLID,3,BLUE);
  pDC->SelectObject(&CutLPen);

  for( i = 0; i < nLine ; i++)
  {
    pVW->GetRE()->DrawLine3DOnClientWindow(pDC,m_arCuttingEdge[i].m_Pos1,m_arCuttingEdge[i].m_Pos2);
  }

  pDC->RestoreDC(-1);
}
*/


