.PHONY: all

all: ${TARGET_DIR}/lesv.sh \
	${TARGET_DIR}/x_hqx2callsv.sh \
	${TARGET_DIR}/x_hqx2makecfg.sh \
	${TARGET_DIR}/x_hqx2splitseq.sh

${TARGET_DIR}/lesv.sh: scripts/lesv.sh
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR}/ ; fi
	cp -pf scripts/lesv.sh ${TARGET_DIR}/lesv.sh
	chmod +x ${TARGET_DIR}/lesv.sh

${TARGET_DIR}/x_hqx2callsv.sh: scripts/x_hqx2callsv.sh
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR}/ ; fi
	cp -pf scripts/x_hqx2callsv.sh ${TARGET_DIR}/x_hqx2callsv.sh
	chmod +x ${TARGET_DIR}/x_hqx2callsv.sh

${TARGET_DIR}/x_hqx2makecfg.sh: scripts/x_hqx2makecfg.sh
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR}/ ; fi
	cp -pf scripts/x_hqx2makecfg.sh ${TARGET_DIR}/x_hqx2makecfg.sh
	chmod +x ${TARGET_DIR}/x_hqx2makecfg.sh

${TARGET_DIR}/x_hqx2splitseq.sh: scripts/x_hqx2splitseq.sh
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR}/ ; fi
	cp -pf scripts/x_hqx2splitseq.sh ${TARGET_DIR}/x_hqx2splitseq.sh
	chmod +x ${TARGET_DIR}/x_hqx2splitseq.sh

clean:
	rm -f ${TARGET_DIR}/lesv.sh
	rm -f ${TARGET_DIR}/x_hqx2callsv.sh
	rm -f ${TARGET_DIR}/x_hqx2makecfg.sh
	rm -f ${TARGET_DIR}/x_hqx2splitseq.sh
