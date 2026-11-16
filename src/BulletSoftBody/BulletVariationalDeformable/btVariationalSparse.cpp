/*
Bullet Continuous Collision Detection and Physics Library
https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "btVariationalSparse.h"

void btVariationalSparseMatrix::setFromTriplets(int n, const btAlignedObjectArray<int>& rows, const btAlignedObjectArray<int>& cols, const btAlignedObjectArray<double>& vals)
{
	m_n = n;
	const int nt = rows.size();
	btAlignedObjectArray<int> start;
	start.resize(n + 1, 0);
	for (int t = 0; t < nt; t++) start[rows[t] + 1]++;
	for (int r = 0; r < n; r++) start[r + 1] += start[r];
	btAlignedObjectArray<int> fill, tc;
	btAlignedObjectArray<double> tv;
	fill.resize(n);
	for (int r = 0; r < n; r++) fill[r] = start[r];
	tc.resize(nt);
	tv.resize(nt);
	for (int t = 0; t < nt; t++)
	{
		int p = fill[rows[t]]++;
		tc[p] = cols[t];
		tv[p] = vals[t];
	}
	m_rowStart.resize(n + 1);
	m_cols.resize(0);
	m_vals.resize(0);
	m_cols.reserve(nt);
	m_vals.reserve(nt);
	for (int r = 0; r < n; r++)
	{
		m_rowStart[r] = m_cols.size();
		const int b = start[r], e = start[r + 1];
		// rows are short, an insertion sort by column is enough
		for (int i = b + 1; i < e; i++)
		{
			int kc = tc[i];
			double kv = tv[i];
			int j = i - 1;
			while (j >= b && tc[j] > kc)
			{
				tc[j + 1] = tc[j];
				tv[j + 1] = tv[j];
				j--;
			}
			tc[j + 1] = kc;
			tv[j + 1] = kv;
		}
		for (int i = b; i < e; i++)
		{
			int last = m_cols.size() - 1;
			if (last >= m_rowStart[r] && m_cols[last] == tc[i])
				m_vals[last] += tv[i];
			else
			{
				m_cols.push_back(tc[i]);
				m_vals.push_back(tv[i]);
			}
		}
	}
	m_rowStart[n] = m_cols.size();
}

void btVariationalSparseMatrix::multiply(const btVec3d* x, btVec3d* y) const
{
	for (int r = 0; r < m_n; r++)
	{
		btVec3d s(0, 0, 0);
		for (int p = m_rowStart[r]; p < m_rowStart[r + 1]; p++) s += x[m_cols[p]] * m_vals[p];
		y[r] = s;
	}
}

namespace
{
///Nested dissection ordering from breadth-first level structures (George's automatic nested dissection)
struct btNestedDissection
{
	const btVariationalSparseMatrix& m_A;
	btAlignedObjectArray<int>& m_order;  // elimination order, filled in
	btAlignedObjectArray<int> m_nodes;   // working permutation of all vertices
	btAlignedObjectArray<int> m_region;  // region stamp per vertex
	btAlignedObjectArray<int> m_visit;   // visit stamp per vertex
	btAlignedObjectArray<int> m_level;
	btAlignedObjectArray<int> m_queue;
	btAlignedObjectArray<int> m_tmp;
	int m_stamp;
	int m_visitStamp;

	btNestedDissection(const btVariationalSparseMatrix& A, btAlignedObjectArray<int>& order) : m_A(A), m_order(order), m_stamp(0), m_visitStamp(0)
	{
		const int n = A.m_n;
		m_nodes.resize(n);
		m_region.resize(n, 0);
		m_visit.resize(n, 0);
		m_level.resize(n, 0);
		m_queue.reserve(n);
		m_tmp.reserve(n);
		for (int i = 0; i < n; i++) m_nodes[i] = i;
		m_order.resize(0);
		m_order.reserve(n);
	}

	///breadth-first search restricted to the current region; fills m_queue (visit order) and m_level
	int bfs(int start, int region)
	{
		m_visitStamp++;
		m_queue.resize(0);
		m_queue.push_back(start);
		m_visit[start] = m_visitStamp;
		m_level[start] = 0;
		int head = 0;
		while (head < m_queue.size())
		{
			int v = m_queue[head++];
			for (int p = m_A.m_rowStart[v]; p < m_A.m_rowStart[v + 1]; p++)
			{
				int w = m_A.m_cols[p];
				if (m_region[w] != region || m_visit[w] == m_visitStamp) continue;
				m_visit[w] = m_visitStamp;
				m_level[w] = m_level[v] + 1;
				m_queue.push_back(w);
			}
		}
		return m_level[m_queue[m_queue.size() - 1]] + 1;  // number of levels
	}

	void emit(int begin, int end)
	{
		for (int i = begin; i < end; i++) m_order.push_back(m_nodes[i]);
	}

	void dissect(int begin, int end)
	{
		const int leafSize = 40;
		while (end - begin > leafSize)
		{
			const int region = ++m_stamp;
			for (int i = begin; i < end; i++) m_region[m_nodes[i]] = region;
			// pseudo-peripheral start vertex: the far end of a first search
			bfs(m_nodes[begin], region);
			int numLevels = bfs(m_queue[m_queue.size() - 1], region);
			const int reached = m_queue.size();
			if (reached < end - begin)
			{
				// the region is disconnected: order this component on its own, continue with the rest
				m_tmp.resize(0);
				for (int i = 0; i < reached; i++) m_tmp.push_back(m_queue[i]);
				for (int i = begin; i < end; i++)
					if (m_visit[m_nodes[i]] != m_visitStamp) m_tmp.push_back(m_nodes[i]);
				for (int i = begin; i < end; i++) m_nodes[i] = m_tmp[i - begin];
				dissect(begin, begin + reached);
				begin += reached;
				continue;
			}
			if (numLevels < 3) break;
			// choose the level that balances the two sides best
			m_tmp.resize(0);
			m_tmp.resize(numLevels, 0);
			for (int i = 0; i < reached; i++) m_tmp[m_level[m_queue[i]]]++;
			int best = 1, below = m_tmp[0];
			int bestImbalance = end - begin;
			for (int l = 1; l < numLevels - 1; l++)
			{
				int above = (end - begin) - below - m_tmp[l];
				int imbalance = below > above ? below - above : above - below;
				if (imbalance < bestImbalance)
				{
					bestImbalance = imbalance;
					best = l;
				}
				below += m_tmp[l];
			}
			// separator: vertices of the chosen level that touch the next level
			int nA = 0, nB = 0, nS = 0;
			m_tmp.resize(0);
			m_tmp.resize(reached, 0);  // 0 = A, 1 = B, 2 = separator, indexed like m_queue
			for (int i = 0; i < reached; i++)
			{
				int v = m_queue[i];
				int side = 0;
				if (m_level[v] > best)
					side = 1;
				else if (m_level[v] == best)
				{
					for (int p = m_A.m_rowStart[v]; p < m_A.m_rowStart[v + 1]; p++)
					{
						int w = m_A.m_cols[p];
						if (m_region[w] == region && m_level[w] == best + 1)
						{
							side = 2;
							break;
						}
					}
				}
				m_tmp[i] = side;
				if (side == 0) nA++;
				if (side == 1) nB++;
				if (side == 2) nS++;
			}
			if (nA == 0 || nB == 0) break;
			int pos[3] = {begin, begin + nA, begin + nA + nB};
			for (int i = 0; i < reached; i++) m_nodes[pos[m_tmp[i]]++] = m_queue[i];
			// the queue is reused by the recursive calls, so remember the split before descending
			const int aEnd = begin + nA, bEnd = begin + nA + nB;
			dissect(begin, aEnd);
			dissect(aEnd, bEnd);
			emit(bEnd, end);
			return;
		}
		emit(begin, end);
	}
};
}  // namespace

void btVariationalSparseCholesky::computeOrdering(const btVariationalSparseMatrix& A)
{
	btNestedDissection nd(A, m_perm);
	nd.dissect(0, A.m_n);
	m_invPerm.resize(A.m_n);
	for (int k = 0; k < A.m_n; k++) m_invPerm[m_perm[k]] = k;
}

bool btVariationalSparseCholesky::factorize(const btVariationalSparseMatrix& A)
{
	const int n = A.m_n;
	m_n = n;
	computeOrdering(A);

	// upper triangle of the permuted matrix, stored by columns
	btAlignedObjectArray<int> cp, ci;
	btAlignedObjectArray<double> cx;
	cp.resize(n + 1);
	ci.reserve(A.nonZeros() / 2 + n);
	cx.reserve(A.nonZeros() / 2 + n);
	for (int k = 0; k < n; k++)
	{
		cp[k] = ci.size();
		const int old = m_perm[k];
		for (int p = A.m_rowStart[old]; p < A.m_rowStart[old + 1]; p++)
		{
			int i = m_invPerm[A.m_cols[p]];
			if (i <= k)
			{
				ci.push_back(i);
				cx.push_back(A.m_vals[p]);
			}
		}
	}
	cp[n] = ci.size();

	// elimination tree
	btAlignedObjectArray<int> parent, ancestor;
	parent.resize(n);
	ancestor.resize(n);
	for (int k = 0; k < n; k++)
	{
		parent[k] = -1;
		ancestor[k] = -1;
		for (int p = cp[k]; p < cp[k + 1]; p++)
		{
			int i = ci[p];
			while (i != -1 && i < k)
			{
				int next = ancestor[i];
				ancestor[i] = k;
				if (next == -1) parent[i] = k;
				i = next;
			}
		}
	}

	// symbolic factorization: count the entries of each column of L
	btAlignedObjectArray<int> flag, next;
	flag.resize(n);
	m_colStart.resize(n + 1);
	next.resize(n);
	for (int k = 0; k < n; k++)
	{
		next[k] = 1;  // the diagonal
		flag[k] = -1;
	}
	for (int k = 0; k < n; k++)
	{
		flag[k] = k;
		for (int p = cp[k]; p < cp[k + 1]; p++)
		{
			for (int i = ci[p]; flag[i] != k; i = parent[i])
			{
				next[i]++;
				flag[i] = k;
			}
		}
	}
	m_colStart[0] = 0;
	for (int k = 0; k < n; k++) m_colStart[k + 1] = m_colStart[k] + next[k];
	m_rowIdx.resize(m_colStart[n]);
	m_Lx.resize(m_colStart[n]);

	// numeric up-looking factorization
	btAlignedObjectArray<int> stack;
	btAlignedObjectArray<double> x;
	stack.resize(n);
	x.resize(n, 0.0);
	for (int k = 0; k < n; k++)
	{
		next[k] = m_colStart[k];
		flag[k] = -1;
	}
	for (int k = 0; k < n; k++)
	{
		// nonzero pattern of row k of L, in topological order, on top of the stack
		int top = n;
		flag[k] = k;
		for (int p = cp[k]; p < cp[k + 1]; p++)
		{
			int i = ci[p];
			x[i] = cx[p];
			int len = 0;
			for (; flag[i] != k; i = parent[i])
			{
				stack[len++] = i;
				flag[i] = k;
			}
			while (len > 0) stack[--top] = stack[--len];
		}
		double d = x[k];
		x[k] = 0.0;
		for (; top < n; top++)
		{
			const int i = stack[top];
			const double lki = x[i] / m_Lx[m_colStart[i]];
			x[i] = 0.0;
			for (int p = m_colStart[i] + 1; p < next[i]; p++) x[m_rowIdx[p]] -= m_Lx[p] * lki;
			d -= lki * lki;
			const int q = next[i]++;
			m_rowIdx[q] = k;
			m_Lx[q] = lki;
		}
		if (!(d > 0.0)) return false;
		const int q = next[k]++;
		m_rowIdx[q] = k;
		m_Lx[q] = sqrt(d);
	}
	m_work.resize(n);
	return true;
}

void btVariationalSparseCholesky::solve(const btVec3d* b, btVec3d* x) const
{
	const int n = m_n;
	btVec3d* y = &m_work[0];
	for (int k = 0; k < n; k++) y[k] = b[m_perm[k]];
	for (int j = 0; j < n; j++)
	{
		const int s = m_colStart[j], e = m_colStart[j + 1];
		y[j] *= 1.0 / m_Lx[s];
		const btVec3d yj = y[j];
		for (int p = s + 1; p < e; p++) y[m_rowIdx[p]] -= yj * m_Lx[p];
	}
	for (int j = n - 1; j >= 0; j--)
	{
		const int s = m_colStart[j], e = m_colStart[j + 1];
		btVec3d yj = y[j];
		for (int p = s + 1; p < e; p++) yj -= y[m_rowIdx[p]] * m_Lx[p];
		y[j] = yj * (1.0 / m_Lx[s]);
	}
	for (int k = 0; k < n; k++) x[m_perm[k]] = y[k];
}
