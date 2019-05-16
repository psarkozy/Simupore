import numpy as np
import time
import os

import parasail

def parasail_cigar_int_to_op(pscigint):
    return ("MIDNSHP=XB"[pscigint & 0xf],  pscigint>>4)



def align_local_maxscore(subject, query, match_score = 1, mismatch_score = -1, deletion_score = -1, insertion_score = -1, onlyscore = False):
    M = np.zeros((len(subject) + 1, len(query) + 1), dtype=np.int32, order= 'C') #int16 means dont go over more than 10k long seq's
    for i in xrange(len(subject) + 1):
        M[i, 0] = 0
    for i in xrange(len(query) + 1):
        M[0, i] = 0
    deletion_score = np.array(deletion_score,dtype = np.int32)
    insertion_score = np.array(insertion_score, dtype = np.int32)
    match_score = np.array(match_score, dtype = np.int32)
    mismatch_score = np.array(mismatch_score, dtype = np.int32)
    zero = np.zeros((1))
    for x, subject_char in enumerate(subject):
        for y, query_char in enumerate( query):
            if subject_char == query_char:
                ms = match_score
            else:
                ms = mismatch_score
            M[x + 1, y + 1] = max(M[x, y] + ms,
                                  M[x, y + 1] + deletion_score,
                                  M[x + 1, y] + insertion_score,
                                  0)
    if onlyscore:
        return np.amax(M)
    else:
        return backtrack_local(M,subject,query)


def align_local_maxscore2(subject, query, match_score = 1, mismatch_score = -1, deletion_score = -1, insertion_score = -1, onlyscore = False):
    M = np.zeros((len(subject) + 1, len(query) + 1), dtype=np.int32, order= 'F') #int16 means dont go over more than 10k long seq's
    for i in xrange(len(subject) + 1):
        M[i, 0] = 0
    for i in xrange(len(query) + 1):
        M[0, i] = 0

    query_ord = np.array(map(ord,query),dtype=np.int32)

    for x, subject_char in enumerate(subject):
        #matchvec = query_ord.where(x == ord(subject_char), match_score,mismatch_score)
        rowmat = np.zeros((4,len(query)+1))

        rowmat[0,:] = 0
        rowmat[:,0] = 0

        rowmat[1,:] = M[x,:] + deletion_score  # ins_score
        rowmat[2,:] = M[x,:] + deletion_score  # delscore
        rowmat[3,:] = M[x,:] + deletion_score  # match score

        M[x+1,:] = np.amax(rowmat,axis = 0)

        for y, query_char in enumerate( query):
            if subject_char == query_char:
                ms = match_score
            else:
                ms = mismatch_score
            M[x + 1, y + 1] = max(M[x, y] + ms,
                                  M[x, y + 1] + deletion_score,
                                  M[x + 1, y] + insertion_score,
                                  0)
    if onlyscore:
        return np.amax(M)
    else:
        return backtrack_local(M,subject,query)


def backtrack_local(M,subject,query):
    bestpos = np.unravel_index(np.argmax(M),M.shape)
    bestscore = np.amax(M)
    currscore = bestscore
    currpos = bestpos
    alignment = {'matches':0, 'mismatches':0,'insertions':0,'deletions':0,'maxscore':bestscore,'bestpos':bestpos}
    #also generate alinged seq's?
    aligned_subject = [subject[currpos[0] -1 ]]
    aligned_query   = [query  [currpos[1] -1 ]]

    while (currscore > 0):
        upleft =  M[currpos[0]-1,currpos[1]-1]
        left =  M[currpos[0]-1,currpos[1]] #deletion
        up =  M[currpos[0],currpos[1]-1] #insertion
        #if max(upleft,left,up) == 0:
        #    break
        if upleft >= max(left,up):
            if  upleft < currscore: #match
                alignment['matches']+=1
                aligned_subject.append(subject[currpos[0]- 1 - 1])
                aligned_query.append(query  [currpos[1] - 1 - 1])
            else:
                alignment['mismatches']+=1
                aligned_subject.append(subject[currpos[0] - 1 - 1].lower())
                aligned_query.append(query  [currpos[1] - 1 - 1].lower())
            currscore = upleft
            currpos = (currpos[0]-1,currpos[1]-1)

        else:
            if left>=up:
                alignment['deletions']+=1
                currscore = up
                currpos = (currpos[0]-1,currpos[1])
                aligned_subject.append(subject[currpos[0] - 1 ])
                aligned_query.append('-')
            else:
                alignment['insertions']+=1
                currscore = up
                currpos = (currpos[0],currpos[1]-1)
                aligned_subject.append('-')
                aligned_query.append(query[currpos[1] - 1])
    alignment['aligned_query'] = ''.join(aligned_query)[::-1] #reverse
    alignment['aligned_subject'] = ''.join(aligned_subject)[::-1] #reverse
    return alignment


class NeedlemanWunsch:
    def __init__(self,substitution_matrix = None,insertion_penalty = -1,deletion_penalty = -1):
        self.substitution = substitution_matrix
        self.insertion_penalty = insertion_penalty
        self.deletion_penalty = deletion_penalty
    def ez_defaults(self):
        self.substitution = {}
        for a in 'ACGT':
            for b in 'ACGT':
                if a==b:
                    self.substitution[(a,b)] = 1
                else:
                    self.substitution[(a,b)] = -1
        self.insertion_penalty = -1
        self.deletion_penalty = -1


    def align_global(self, subject,query):
        self.subject = subject
        self.query = query
        dp = self.deletion_penalty
        ip = self.insertion_penalty

        M = np.zeros((len(subject) + 1, len(query) + 1))
        for i in xrange(len(subject) + 1):
            M[i, 0] = -i
        for i in xrange(len(query) + 1):
            M[0, i] = -i
        #for x in xrange(1, len(subject)+1):
            #print x, len(subject)
            #for y in xrange(1, len(query) +1):
                #M[x,y] = max(M[x-1,y-1] + self.substitution[(subject[x-1],query[y-1])], M[x-1,y] +self.deletion_penalty, M[x,y-1] + self.insertion_penalty )
            #    M[x,y] = max(M[x-1,y-1] + (1 if subject[x-1] == query[y-1] else -1), M[x-1,y] +self.deletion_penalty, M[x,y-1] + self.insertion_penalty )
        for x, subject_char in enumerate(subject):
            for y, query_char in enumerate(query):
                M[x + 1, y + 1] = max(M[x, y] + (1 if subject_char == query_char else -1), M[x,y+1] + dp, M[x+1,y] + ip)

        self.M=M

    def align_local(self, subject,
                    query):  # only halfassedly implemented, thus only returns max score (the rest doesnt really matter for smith waterman anyway)
        self.subject = subject
        self.query = query

        dp = self.deletion_penalty
        ip = self.insertion_penalty

        M = np.zeros((len(subject) + 1, len(query) + 1), dtype=np.int32)
        for i in xrange(len(subject) + 1):
            M[i, 0] = 0
        for i in xrange(len(query) + 1):
            M[0, i] = 0
        # for x in xrange(1, len(subject)+1):
        #    for y in xrange(1, len(query) +1):
        #        M[x,y] = max(M[x-1,y-1] + self.substitution[(subject[x-1],query[y-1])], M[x-1,y] +self.deletion_penalty, M[x,y-1] + self.insertion_penalty,0 )

        for x, subject_char in enumerate(subject):
            for y, query_char in enumerate(query):
                M[x + 1, y + 1] = max(M[x, y] + (1 if subject_char == query_char else -1), M[x, y + 1] + dp,
                                      M[x + 1, y] + ip, 0)

        for x in xrange(1, len(subject) + 1):
            for y in xrange(1, len(query) + 1):
                M[x, y] = max(M[x - 1, y - 1] + self.substitution[(subject[x - 1], query[y - 1])],
                              M[x - 1, y] + self.deletion_penalty, M[x, y - 1] + self.insertion_penalty, 0)
        self.M = M
        return np.amax(M)


    def backtrack_global(self):
        insertions = 0
        deletions = 0
        matches = 0
        mismatches = 0
        score = 0
        current_x = self.M.shape[0] - 1
        current_y = self.M.shape[1] - 1
        qstr = ''
        sstr = ''
        while current_x > 0 and current_y > 0:
            if self.M[current_x-1,current_y-1] >= self.M[current_x, current_y -1] and self.M[current_x-1,current_y-1] >= self.M[current_x - 1, current_y]:
                if self.subject[current_x-1] == self.query[current_y-1]:
                    matches += 1
                    sstr = self.subject[current_x - 1] + sstr
                    qstr = self.query[current_y - 1] +   qstr
                else:
                    mismatches +=1
                    sstr = str.lower(self.subject[current_x - 1]) + sstr
                    qstr = str.lower(self.query[current_y - 1]) +   qstr
                score += self.substitution[(self.subject[current_x-1], self.query[current_y-1])]
                current_y -= 1
                current_x -= 1
            elif self.M[current_x, current_y -1] <= self.M[current_x -1, current_y]: #deletion
                qstr = '-' + qstr
                sstr = self.subject[current_x - 1] + sstr
                deletions += 1
                score += self.deletion_penalty
                current_x -= 1
            else:

                qstr = self.query[current_y - 1] + qstr
                sstr = '-'+ sstr
                insertions += 1
                score += self.insertion_penalty
                current_y -= 1

        score += self.insertion_penalty * current_y
        score += self.deletion_penalty *  current_x
        deletions += current_x
        insertions += current_y
        print sstr
        print qstr
        return score,matches,mismatches,insertions,deletions



    def backtrack_local(self):
        insertions = 0
        deletions = 0
        matches = 0
        mismatches = 0
        score = 0
        bestalignstart= np.argmax(self.M)
        current_x = bestalignstart[0]
        current_y = bestalignstart[1]

        while score >0:
            if self.M[current_x-1,current_y-1] >= self.M[current_x, current_y -1] and self.M[current_x-1,current_y-1] >= self.M[current_x - 1, current_y]:
                if self.subject[current_x] == self.query[current_y]:
                    matches += 1
                else:
                    mismatches +=1
                score += self.substitution[(self.subject[current_x], self.query[current_y])]
                current_y -= 1
                current_x -= 1
            elif self.M[current_x, current_y -1] <= self.M[current_x -1, current_y]: #deletion
                deletions += 1
                score += self.deletion_penalty
                current_x -= 1
            else:
                insertions += 1
                score += self.insertion_penalty
                current_y -= 1

        score += self.insertion_penalty * current_y
        score += self.deletion_penalty *  current_x
        deletions += current_x
        insertions += current_y

        return score,matches,mismatches,insertions,deletions

def tracefromcigarlist(query, subject, cigars):
    str1 = ''
    str2 = ''
    cigarlength = 0
    for (cigop, l) in cigars:
        # print 'cigarlength = %d op = %s Len(str1)=%d, len(str2)=%d len(query)=%d len(subject)=%d totalops = %d'%(l,cigop,len(str1),len(str2),len(query),len(subject),cigarlength)
        cigarlength += l
        if cigop == '=':
            str1 += subject[0:l]
            subject = subject[l:]
            str2 += query[0:l]
            query = query[l:]
        elif cigop == 'X':
            if subject[0:l] != query[0:l]:
                str1 += subject[0:l].lower()
                str2 += query[0:l].lower()
            else:
                str1 += subject[0:l]
                str2 += query[0:l]

            subject = subject[l:]
            query = query[l:]
        elif cigop == 'D':
            str1 += subject[0:l]
            subject = subject[l:]
            str2 += '-' * l
            # query = query[l:]
        elif cigop == 'I':
            str1 += '+' * l
            # subject = subject[l:]
            str2 += query[0:l]
            query = query[l:]

    for i in range(0, len(str1), 80):
        print i, str1[i:min(len(str1), i + 80)]
        print i, str2[i:min(len(str2), i + 80)]
    print 'cigarlength', cigarlength, 'strlens:', len(str1), len(str2)
    return str1,str2


def parasail_nw_onlyscore(query, subject, printtrace=False,verbose =False):
    t0 = time.clock()
    mat = parasail.matrix_create('ACGT',1,-1)
    result = parasail.nw_trace_diag_32(query,subject,1,1,mat)
    #result2 = parasail.sg_trace_scan(query,subject,1,1,mat)
    cigseq = result.cigar.seq
    cigarops = map(parasail_cigar_int_to_op,cigseq)
    qpos = 0
    spos = 0
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    for (op,l) in cigarops:
        if op == '=':
            qpos += l
            spos += l
            matches += l
        elif op == 'D':
            spos += l
            deletions += l
        elif op == 'I':
            qpos += l
            insertions += l
        elif op == 'X':
            for c in xrange(l):
                if query[qpos + c] == subject[spos +c] :
                    matches += 1
                else:
                    mismatches +=1
            qpos += l
            spos += l



    if printtrace:
        tracefromcigarlist(query,subject,cigarops)
    #print result2.score
    #print result2.cigar.seq
    #print map(parasail_cigar_int_to_op, result2.cigar.seq)
 #cig = result.cigar
    t0 = time.clock()-t0
    output =  {'score':result.score,'matches:':matches,'mismatches':mismatches,'insertions':insertions,'deletions':deletions}
    if verbose:
        print 'ParasailNW Q%d to S%d score is %d t=%.3f'%(len(query),len(subject),result.score,t0)
        print output
    return output




if __name__ == "__main__":
    s = 'GGGGGGG'+('ATTACAGATTACA' * 200)
    q = ('CATTAGAAACATTAGA' * 200) + ''
    #s = 'ACGTACGTACGT'
    #q = 'ACGTCGGACGT'
    s = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACATTAGAAACATTAGACATTAGAAACATTAGACATTAGAAACATTAGAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    q = 'CATTAGAAACATTAGACATTAGAAACATTAGACATTAGAAACATTAGA'

    t0 = time.clock()
    mat = parasail.matrix_create('ACGT',1,-1)
    #result = parasail.sw_trace_striped_16(s,q,1,1,mat)
    result = parasail.sw_trace_diag_32(q,s,1,1,mat)
    cig = result.cigar
#    print cig.seq
    print(result.end_ref)
    print(result)
    #print(result.read_begin1)
    #print(result.read_end1)
    print map(parasail_cigar_int_to_op, cig.seq)

    #print cig.decode()

 #   print( result.cigar)
    # nw.align_global(s,q)
    print 'Time taken for parasail',len(s),'x',len(q),':',time.clock() -t0,'score = ',result.score


    t0 = time.clock()
    score = align_local_maxscore(s, q,onlyscore=True)
    # nw.align_global(s,q)
    print 'Time taken for',len(s),'x',len(q),':',time.clock() -t0,'score = ',score

    #exit(1)
    t0 = time.clock()
    score = align_local_maxscore2(s, q,onlyscore=True)
    # nw.align_global(s,q)
    print 'Time taken for',len(s),'x',len(q),':',time.clock() -t0,'score = ',score

    exit(1)

    nw = NeedlemanWunsch()
    nw.ez_defaults()

   # print 'Align time = %f, Backtrack time = %f'%(t1-t0,time.clock()-t1)
    print nw.backtrack_global()
    t0=time.clock()
    os.system('c:/anaconda2/scripts/nwalign.exe %s %s'%(s+s,q+q))
    print time.clock()-t0

