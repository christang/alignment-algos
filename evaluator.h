/**
 *  Package HMAP2.1
 *  File: evaluator.h
 *  Desc: Base template class for defining evaluators.
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_EVALUATOR
#define _HMAP2_EVALUATOR

#include <string>

class SimilarityMatrix;

template <class S1, class S2, class Etype>
class Evaluator {

public:
  
  /**
   *  Similarity function for an alignment.  To be defined by derived class!
   *  Defines the similarity between one element in the query and
   *  one in the template sequence.  To be defined by derived class!
   * @return Specifies similarity of position, to be maximized by alignment.
   */
  float similarity (const S1& query_seq,
		    const S2& templ_seq,
		    int position_in_query,
		    int position_in_templ
		    ) const ;
  
  /**
   *  Deletion penalty for an alignment.  To be defined by derived class!
   *  Defines the deletion penalty for removing sequence elements
   *  between positions t1 to t2 in the template from the query 
   *  while aligning q1 to t1 and q2 to t2.  To be defined by derived class!
   * @param position1_in_query aligns to position1_in_templ, pre-deletion
   * @param position2_in_query aligns to position2_in_templ, post-deletion
   * @param position1_in_templ aligns to position1_in_query
   * @param position2_in_templ aligns to position2_in_query
   * @return Specifies penalty to be minimized by alignment.
   */
  float   deletion (const S1& query_seq,
		    const S2& templ_seq,
		    int position1_in_query,
		    int position2_in_query,
		    int position1_in_templ,
		    int position2_in_templ
		    ) const ;
  /**
   *  Insertion penalty for an alignment.  To be defined by derived class!
   *  Defines the insertion penalty for having extra sequence
   *  elements between positions q1 to q2 in the query 
   *  while aligning q1 to t1 and q2 to t2.  To be defined by derived class!
   * @param position1_in_query aligns to position1_in_templ, pre-insertion
   * @param position2_in_query aligns to position2_in_templ, post-insertion
   * @param position1_in_templ aligns to position1_in_query
   * @param position2_in_templ aligns to position2_in_query
   * @return Specifies penalty to be minimized by alignment.
   */
  float  insertion (const S1& query_seq,
		    const S2& templ_seq,
		    int position1_in_query,
		    int position2_in_query,
		    int position1_in_templ,
		    int position2_in_templ
		    ) const ;
  
  void post_process (SimilarityMatrix& s) const ;
  void pre_calculate (const S1& query_seq, const S2& templ_seq) const;

protected:
  
  /**
   *  Constructor is protected.
   *  Evaluators are defined only within derived classes and
   *  only derived classes may create instances of an evaluator.
   */
  Evaluator ()  {};

  /**
   *  Returns '*this' recast as the derived class.
   *  This trick avoids the use of virtual functions to define
   *  the behavior of similarity(), insertion() and deletion() in
   *  the derived class, which is execution inefficient.  Instead,
   *  the strategy is to use the Barton & Nackman trick as 
   *  discussed here: Curiously recursive template pattern, 
   *  http://osl.iu.edu/~tveldhui/papers/techniques/
   */
  inline const Etype& Derived() const { return static_cast<const Etype&>(*this); } ;
  
};


template <class S1, class S2, class Etype>
inline float Evaluator<S1,S2,Etype>::similarity (const S1& query_seq,
						 const S2& templ_seq,
						 int position_in_query,
						 int position_in_templ
						 ) const
{ return Derived().similarity(query_seq,templ_seq,
			      position_in_query,position_in_templ); }

template <class S1, class S2, class Etype>
inline float Evaluator<S1,S2,Etype>::deletion (const S1& query_seq,
					       const S2& templ_seq,
					       int position1_in_query,
					       int position2_in_query,
					       int position1_in_templ,
					       int position2_in_templ
					       ) const
{ return Derived().deletion(query_seq,templ_seq,
			    position1_in_query,position2_in_query,
			    position1_in_templ,position2_in_templ); }


template <class S1, class S2, class Etype>
inline float Evaluator<S1,S2,Etype>::insertion (const S1& query_seq,
						const S2& templ_seq,
						int position1_in_query,
						int position2_in_query,
						int position1_in_templ,
						int position2_in_templ
						) const
{ 
  return Derived().insertion(query_seq,templ_seq,
			     position1_in_query,position2_in_query,
			     position1_in_templ,position2_in_templ); }


template <class S1, class S2, class Etype>
inline void Evaluator<S1,S2,Etype>::post_process (SimilarityMatrix& s) const
{
  return Derived().post_process(s);
}

template <class S1, class S2, class Etype>
inline void Evaluator<S1,S2,Etype>::pre_calculate (const S1& query_seq,
                                                   const S2& templ_seq) const
{
  return Derived().pre_calculate(query_seq,templ_seq);
}

						  
#endif  // _HMAP2_EVALUATOR
