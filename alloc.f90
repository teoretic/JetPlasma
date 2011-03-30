module alloc

use names

implicit none

contains

subroutine strt
  integer :: err

  print *
  print *, 'alloc()'
  print *, '==========='

  allocate(U(5,nElems+nElems_b), stat = err)
  print *, 'U allocate stat=', err
  allocate(U_hat(5,nElems+nElems_b), stat = err)
  print *, 'U_hat allocate stat=', err
  allocate(U_til(5,nElems+nElems_b), stat = err)
  print *, 'U_til allocate stat=', err
  allocate(B(3,nElems+nElems_b), stat = err)
  print *, 'B allocate stat=', err
  allocate(B_hat(3,nElems+nElems_b), stat = err)
  print *, 'B_hat allocate stat=', err
  allocate(B_angles(2,3,nElems+nElems_b), stat = err)
  print *, 'B_angles allocate stat=', err
  allocate(B_hat_angles(2,3,nElems+nElems_b), stat = err)
  print *, 'B_hat_angles allocate stat=', err
  allocate(U_angles(2,5,nElems+nElems_b), stat = err)
  print *, 'U_angles allocate stat=', err
  allocate(U_hat_angles(2,5,nElems+nElems_b), stat = err)
  print *, 'U_hat_angles allocate stat=', err
  allocate(gas_edge_flux(1:nEdges,5), stat = err)
  print *, 'gas_edge_flux allocate stat=', err
  allocate(up_koef(2,3,nElems+nElems_b), stat = err)
  print *, 'up_koef allocate stat=', err
  allocate(upstr_nodes(10,nNodes), stat = err)
  print *, 'upstr_nodes allocate stat=', err
  allocate(up_koef_z(3,nElems+nElems_b), stat = err) ! три грани ячейки
  print *, 'up_koef_z allocate stat=', err

  allocate(B_edges(2,nEdges), stat = err)
  print *, 'B_edges allocate stat=', err
  allocate(B_edges_hat(3,nEdges), stat = err)
  print *, 'B_edges_hat allocate stat=', err
  allocate(V_edges(2,nEdges), stat = err)
  print *, 'V_edges allocate stat=', err
  allocate(B_nodes(3,nNodes), stat = err)
  print *, 'B_nodes allocate stat=', err
  allocate(V_nodes(2,nNodes), stat = err)
  print *, 'V_nodes allocate stat=', err
  allocate(gradBB(3,nElems+nElems_b), stat = err)
  print *, 'gradBB allocate stat=', err
  allocate(MagnEn(nElems+nElems_b), stat = err)
  print *, 'MagnEn allocate stat=', err
  allocate(divBB(3,nElems+nElems_b), stat = err)
  print *, 'divBB allocate stat=', err

  allocate(rotB(3,nElems+nElems_b), stat = err)
  print *, 'rotB allocate stat=', err
  allocate(rotB_nodes(nNodes), stat = err)
  print *, 'rotB_nodes allocate stat=', err
  allocate(rotB_edges(nEdges), stat = err)
  print *, 'rotB_edges allocate stat=', err

  print *
  print *, 'alloc() DONE'
  print *, '================'
end subroutine strt

subroutine stp
  deallocate(U,U_hat,U_til,B,B_hat,B_angles,B_hat_angles,U_angles,U_hat_angles,B_edges,B_edges_hat,V_edges,B_nodes,V_nodes,gradBB,MagnEn,divBB,gas_edge_flux)
end subroutine stp

end module alloc