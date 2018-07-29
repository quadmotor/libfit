

$(TORUS_DIR)/main: $(TORUS_DIR)/main.f90 $(LIB_DIR)/libFit.a
	$(FC) $(FFLAGS) $^ -o $@

$(TORUS_DIR)/clean:
	-rm $(TORUS_DIR)/output/* 2>/dev/null || true
